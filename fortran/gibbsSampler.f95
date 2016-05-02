! gibbsSampler.f95
MODULE gibbs_sampler

CONTAINS
      subroutine gibbsSampler(matrix,NZW,NZM,NZ,NM,ntopics, &
     max_iter,M,N,p_z)
        IMPLICIT NONE
! Everything is transpose in fortran (??) so rows and columns get switched
        integer :: M, N,i,j,ll,nn,ntapp
        integer max_iter
        integer Z
        integer matrix(M,N)
        integer ntopics
        real*4 p_z(ntopics)
        integer NZ(ntopics)
        integer NM(M)
        integer NZM(ntopics,M)
        integer NZW(N,ntopics)
        real*4 :: alpha, beta
        EXTERNAL genmul
      ntapp = ntopics-1
      do i=1,max_iter
        do j=1,M
! FIXME: ll is like w in the enumerate
!   (i.e. if first nonzero is at place 27, rep 27 (ll) 4 times (ll)
          do ll=1,N
            if (matrix(j,ll) == 0) then
              cycle
            endif
            ! nn is like i in the enumerate
            do nn=1,matrix(j,ll)
              write(*,*) p_z(1:ntapp) 
              call genmul(1,p_z(1:ntapp),ntopics,Z)
              write(*,*) ll
              write(*,*) nn
              write(*,*) NZM(Z,j)
              NZM(j,Z) = NZM(j,Z) - 1
              NM(j) = NM(j) - 1
              NZW(matrix(j,ll),Z) = NZW(matrix(j,ll),Z) - 1
              NZ(Z) = NZ(Z) - 1

              write(*,*) p_z(1:ntapp)
              write(*,*) sum(p_z)
              call  conditional_distribution(NZW, NZM, NZ, NM, beta, alpha,ntopics,M,N,p_z,j,ll)
              call genmul(1,abs(p_z(1:ntapp)),ntopics,Z)

              NZM(j,Z) = NZM(j,Z) + 1
              NM(j) = NM(j) + 1
              NZW(Z,matrix(j,ll)) = NZW(Z,matrix(ll,z)) + 1
              NZ(Z) = NZ(Z) + 1
            enddo  
          enddo
        enddo
      enddo
      end

! End gibbsSampler



! loglikelihood !/! FIXME: N here needs to be vocab

      subroutine loglikelihood(NZW, NZM, alpha, beta, ntopics,N,M,lik)
        IMPLICIT NONE
        real*4 NZM(ntopics,M)
        real*4 NZW(M,ntopics)
        real*4  :: alpha, beta
        integer :: N, M, z, mm, ntopics
        real*4 lik

        lik = 0d0

        do z = 1,ntopics
              call log_multinomial_beta(NZW(z,:) + beta,0,lik)
              call log_multinomial_beta((/beta/),N,lik)
        enddo

        do mm = 1,ntopics
              call log_multinomial_beta(NZM(mm,:) + alpha,0,lik)
              call log_multinomial_beta((/alpha/),ntopics,lik)
        enddo
      end  

! End loglikelihood      

! log_multinomial_beta
      subroutine log_multinomial_beta(alpha,K,lik)
        IMPLICIT NONE
        real*4, dimension(:) :: alpha
        real*4 lik
        integer K
 
        if (K == 0) then
           lik = lik + sum(log_gamma(alpha)) - log_gamma(sum(alpha))
        else
! in the original code it was -=, but K!=0 only when we subtract
           lik = lik - (K * sum(log_gamma(alpha)) - log_gamma(K * sum(alpha)))
        endif
      end
! End multinomial beta

! Conditional Distribution

      subroutine conditional_distribution(NZW, NZM, NZ, NM, beta, alpha,ntopics,M,N,p_z,j,ll)
      IMPLICIT NONE
        integer :: j,ll, M, N
        integer ntopics
        real*4 p_z(ntopics)
        real*4 beta
        real*4 alpha
        integer NZ(ntopics)
        integer NM(M)
        integer NZM(ntopics,M)
        integer NZW(N,ntopics)
        real*4 left(ntopics)
        real*4 right(ntopics)

        left = (NZW(ll,:) + beta) / (NZ + beta * N)
        right = (NZM(:,j) + alpha) / (NM(j) + alpha * ntopics)
        p_z = abs(left * right)
        p_z = p_z / sum(p_z)
      end
! End Conditional Distribution
END MODULE gibbs_sampler

