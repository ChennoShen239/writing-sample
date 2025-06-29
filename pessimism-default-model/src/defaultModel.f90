!
! Sovereign default model with long-term debt
!
PROGRAM defaultModel
   USE iso_Fortran_env, ONLY: wp => real64
   USE defMod
   IMPLICIT NONE

   INTEGER :: iter
   REAL(wp) :: errV, errQ
   REAL(wp), DIMENSION(bSz) :: cc
   INTEGER :: yIx, bIx, yPrIx, bPrIx

   REAL(wp), DIMENSION(ySz) :: hy, uhy

   REAL(wp) :: Wbar, theSum
   REAL(wp), DIMENSION(bSz) :: theExps
   ! 新增：用于OpenMP并行债券定价的临时变量
   REAL(wp) :: p_lender_default, rhoD_lender, Wbar_lender, theSum_lender
   
   ! 环境变量读取相关的临时变量
   CHARACTER(LEN=100) :: thetaD_str, outDir_env
   INTEGER :: stat

   ! 通过环境变量或命令行参数赋值
   CALL get_environment_variable("THETAD", thetaD_str, status=stat)
   IF (stat == 0 .AND. LEN_TRIM(thetaD_str) > 0) THEN
      READ(thetaD_str, *) thetaD
   ELSE
      thetaD = 1.0_wp  ! 默认baseline
   END IF

   CALL get_environment_variable("OUTDIR", outDir_env, status=stat)
   IF (stat == 0 .AND. LEN_TRIM(outDir_env) > 0) THEN
      outDir = TRIM(outDir_env)
   ELSE
      outDir = "./results/"
   END IF

   WRITE (*, *) "Start! thetaD=", thetaD, " outDir=", TRIM(outDir)
   CALL allocateAll()
   CALL prepareShocksAndGrids()

   hy = hFun(yGrid)
   uhy = uFun(hy)

   ! Initial values
   Vd0 = uhy
   q0 = 1.0_wp

   !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bIx)
   DO yIx = 1, ySz
   DO bIx = 1, bSz
      V0(yIx, bIx) = uFun(MAX(yGrid(yIx) - kappa*bGrid(bIx), 0.01_wp))
   END DO
   END DO

   ! Main lopp
   iter = 1
   errV = 1.0_wp
   errQ = 1.0_wp
   DO WHILE (iter <= maxIter .AND. (errV > tolErrV .OR. errQ > tolErrQ))
      ! V0, Vd0 => Vd1
      !
      !$OMP PARALLEL DO PRIVATE(yIx)
      DO yIx = 1, ySz
         Vd1(yIx) = uhy(yIx) + beta*DOT_PRODUCT(yPi(yIx, :), &
                                                gamm*V0(:, 1) + (1.0_wp - gamm)*Vd0)
      END DO

      ! V0, q0 => W, Vr, bPol
      ! Vd1, Vr => V1, dPol
      !
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bIx,cc,Wbar,theExps,theSum)
      DO yIx = 1, ySz
      DO bIx = 1, bSz
         cc = yGrid(yIx) - kappa*bGrid(bIx) &
              + q0(yIx, :)*(bGrid - (1.0_wp - delta)*bGrid(bIx))

         WHERE (cc <= 0.0_wp)
            W(yIx, bIx, :) = veryNegative
         ELSEWHERE
            W(yIx, bIx, :) = uFun(cc) + beta*MATMUL(yPi(yIx, :), V0)
         END WHERE

         Wbar = MAXVAL(W(yIx, bIx, :))
         theExps = EXP((W(yIx, bIx, :) - Wbar)/rhoB)
         theSum = SUM(theExps)

         Vr(yIx, bIx) = Wbar + rhoB*LOG(theSum)
         bPol(yIx, bIx, :) = theExps/theSum

         Wbar = MAX(Vr(yIx, bIx), Vd1(yIx))
         theSum = EXP((Vd1(yIx) - Wbar)/rhoD) + EXP((Vr(yIx, bIx) - Wbar)/rhoD)
         dPol(yIx, bIx) = EXP((Vd1(yIx) - Wbar)/rhoD)/theSum
         V1(yIx, bIx) = Wbar + rhoD*LOG(theSum)
      END DO
      END DO

      ! q0, dPol, bPol => q1
      !
      !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx, bPrIx, yPrIx, p_lender_default, rhoD_lender, Wbar_lender, theSum_lender)
      DO yIx = 1, ySz
      DO bPrIx = 1, bSz
         q1(yIx, bPrIx) = 0.0_wp
         DO yPrIx = 1, ySz
            ! 1. 定义贷款人"感知"到的冲击参数 (比真实的更大)
            rhoD_lender = thetaD * rhoD

            ! 2. 计算贷款人"认为"的政府违约概率
            !    这个计算使用了与计算dPol完全相同的逻辑，但使用了被夸大的冲击参数rhoD_lender
            Wbar_lender = MAX(Vr(yPrIx, bPrIx), Vd1(yPrIx))
            theSum_lender = EXP((Vd1(yPrIx) - Wbar_lender) / rhoD_lender) + &
                            EXP((Vr(yPrIx, bPrIx) - Wbar_lender) / rhoD_lender)
            
            ! 防止 a/0 的情况 (当 a 和 b 都极小导致 a+b=0)
            IF (theSum_lender > 0.0_wp) THEN
                p_lender_default = EXP((Vd1(yPrIx) - Wbar_lender) / rhoD_lender) / theSum_lender
            ELSE
                p_lender_default = 0.5_wp ! 或者根据情况设为0或1
            END IF

            ! 3. 使用贷款人"感知"的违约概率来更新债券价格
            q1(yIx, bPrIx) = q1(yIx, bPrIx) &
                             + yPi(yIx, yPrIx) * (1.0_wp - p_lender_default) &
                             * (kappa + (1.0_wp - delta) * &
                                DOT_PRODUCT(bPol(yPrIx, bPrIx, :), q0(yPrIx, :)))
         END DO
      END DO
      END DO
      q1 = q1 / (1.0_wp + rf)

      ! Check for convergence and iterate
      errV = MAX( &
             MAXVAL(ABS(V1 - V0)), &
             MAXVAL(ABS(Vd1 - Vd0)))
      errQ = MAXVAL(ABS(q1 - q0))
      IF (MOD(iter, 10) == 0) WRITE (*, "(I5,2ES20.5)") iter, errV, errQ
      iter = iter + 1
      V0 = V1
      Vd0 = Vd1
      q0 = q1
   END DO ! end of main loop WHILE

   !$OMP PARALLEL DO COLLAPSE(2) PRIVATE(yIx,bIx)
   DO yIx = 1, ySz
   DO bIx = 1, bSz
      EbPr(yIx, bIx) = DOT_PRODUCT(bGrid, bPol(yIx, bIx, :))
   END DO
   END DO

   WRITE (*, *) "Saving to disk..."
   CALL saveResults()

   WRITE (*, *) "Simulate..."
   CALL simulate()

   WRITE (*, *) "The end! ^_^"

CONTAINS

END PROGRAM defaultModel
