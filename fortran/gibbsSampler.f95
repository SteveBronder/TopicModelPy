! gibbsSampler.f95
MODULE cool_guy_22

CONTAINS
      subroutine gibbsSampler(matrix,NWZ,NZM,NZ,NM,Z,ntopics, &
     max_iter,M,N,logPw_z)
        IMPLICIT NONE
        integer :: M, N,i,k
        integer max_iter
        integer Z(M,N)
        integer matrix(N,M)
        integer ntopics
        real*8 p(ntopics)
        real*8 NZ(N)
        real*8 NM(M)
        real*8 NZM(ntopics,M)
        real*8 NWZ(M,ntopics)
        real*8 logPw_z(max_iter)
        EXTERNAL genmul
! Everything but comments must be between columns 6 and 76
      do i=1,max_iter
        do m=1,M
! FIXME: This N should probably be document vocabulary not actual words,
!  but the non-zero entries of the doc m
          do n=1,N
            NZM(Z(m,n),m) = NZM(Z(m,n),m) - 1
            NM(m) = NM(m) - 1
            NWZ(matrix(m,n),Z(m,n)) = NWZ(matrix(m,n),Z(m,n)) - 1
            NZ(Z(m,n)) = NZ(Z(m,n)) - 1

            do k=1,ntopics
              p(k) = NWZ(matrix(m,n),k)/NZ(k) * NZM(k,m)
            enddo

            p = p/sum(p)
            call genmul(1,p,ntopics,Z(m,n))

            NZM(Z(m,n),m) = NZM(Z(m,n),m) + 1
            NM(m) = NM(m) + 1
            NWZ(matrix(m,n),Z(m,n)) = NWZ(matrix(m,n),Z(m,n)) + 1
            NZ(Z(m,n)) = NZ(Z(m,n)) + 1
          enddo
        enddo
      enddo
      end

! End gibbsSampler



! loglikelihood !/! FIXME: N here needs to be vocab

      subroutine loglikelihood(NWZ, NZM, alpha, beta, ntopics,N,M,lik)
        IMPLICIT NONE
        real*8 NZM(ntopics,M)
        real*8 NWZ(M,ntopics)
        real*8  :: alpha, beta
        integer :: N, M, z, mm, ntopics
        real*8 lik

        lik = 0d0

        do z = 1,ntopics
              call log_multinomial_beta(NWZ(z,:) + beta,0,lik)
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
        real*8, dimension(:) :: alpha
        real*8 lik
        integer K
! 
!        if (K == 0) then
!           lik = lik + sum(dlgama(alpha)) - dlgama(sum(alpha))
!        else
! in the original code it was -=, but K!=0 only when we subtract
!           lik = lik - (K * sum(dlgama(alpha)) - dlgama(K * sum(alpha)))
!        endif
      end
! End multinomial beta


END MODULE cool_guy_22

