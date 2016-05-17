! gibbsSampler.f95
MODULE gibbs_sampler

CONTAINS
      subroutine gibbsSampler(matrix,NZW,NZM,NZ,NM,ntopics, &
     max_iter,M,N,p_z,topics,topics2, alpha, beta,lik)
        IMPLICIT NONE
! Everything must be row contiguous in fortran
! you can check this with <array>.flags.f_contiguous
        integer*8, dimension(n,m) :: matrix
        integer*8 m
        integer*8 n
        integer*8, dimension(n,ntopics) :: nzw
        integer*8, dimension(ntopics,m) :: nzm
        integer*8, dimension(ntopics) :: nz
        integer*8, dimension(m) :: nm
        integer*8, intent(in) :: ntopics
        integer*4 Genntopics
        integer*8, intent(in) :: max_iter
        integer*8, dimension(n,m) :: topics
        integer*8, dimension(n,m) :: topics2
        integer*8 :: i,j,ll,nn,ntapp
        integer*8 Z
        integer*4 genZ
        real*8, dimension(ntopics) :: p_z
        real*4, dimension(ntopics) :: genp_z
        real*8,intent(in) :: alpha
        real*8,intent(in) :: beta
        real*8, dimension(max_iter), intent(inout) :: lik
        EXTERNAL genmul
      ntapp = ntopics-1
      do i=1,max_iter
        do j=1,M
          do ll=1,N
              Z = topics(ll,j)
              NZM(Z,j) = NZM(Z,j) - 1
              NM(j) = NM(j) - 1
              NZW(n,Z) = NZW(n,Z) - 1
              NZ(Z) = NZ(Z) - 1

              call  conditional_distribution(matrix,NZW, NZM, NZ, beta, alpha,ntopics,M,N,p_z,j,ll)
              
              ! genmul only accebts kind = 4 variables
              genp_z = real(p_z,4)
              genntopics = int(ntopics,4)
              genz = int(Z,4)

              call genmul(1,abs(genp_z(1:ntapp)),genntopics,genZ)
              genZ = genZ + 1
              topics2(ll,j) = genZ
              NZM(Z,j) = NZM(genZ,j) + 1
              NM(j) = NM(j) + 1
              NZW(n,Z) = NZW(n,genZ) + 1
              NZ(Z) = NZ(genZ) + 1
            enddo  
          enddo
        enddo
        write(*,*) 'Iteration:'
        write(*,*) i
        call loglikelihood(matrix,NZW, NZM, alpha, beta, ntopics,N,M,lik,max_iter,i,j)
        write(*,*) 'Likelihood:'
        write(*,*) lik(i)
      enddo
      end

! End gibbsSampler



! loglikelihood !/! FIXME: N here needs to be vocab

      subroutine loglikelihood(matrix,NZW, NZM, alpha, beta, ntopics,N,M,lik,max_iter,i,j)
        IMPLICIT NONE
        integer*8,dimension(n,m),intent(in) :: matrix
        integer*8 vsize
        integer*8,intent(in) :: N, M, ntopics,max_iter,i,j
        integer*8 :: nn,z,mm
        integer*8,dimension(ntopics,M),intent(in) :: NZM
        integer*8,dimension(N,ntopics),intent(in) :: NZW
        real*8, intent(in)  :: alpha, beta
        real*8,dimension(max_iter), intent(inout) :: lik
        vsize = 0
        do nn = 1,N
          if (matrix(nn,j) == 0) then
            cycle
          endif
          vsize = vsize + 1
        enddo

        do z = 1,ntopics
              call log_multinomial_beta(NZW(z,:) + beta,lik,ntopics,max_iter,i)
              ! Because below only takes in a single value we write seperate function
              call log_multinomial_beta_single(beta,vsize,lik,max_iter,i)

        enddo

        do mm = 1,M
              call log_multinomial_beta(NZM(mm,:) + alpha,lik,ntopics,max_iter,i)
              call log_multinomial_beta_single(alpha,ntopics,lik,max_iter,i)
        enddo

      end  

! End loglikelihood      

! log_multinomial_beta
      subroutine log_multinomial_beta(alpha,lik,ntopics,max_iter,i)
        IMPLICIT NONE
        real*8,dimension(ntopics), intent(in) :: alpha
        real*8,dimension(max_iter),intent(inout) :: lik
        integer*8,intent(in) :: ntopics,max_iter,i

           lik(i) = lik(i) + sum(log_gamma(alpha)) - log_gamma(sum(alpha))
      end
     
! End multinomial beta

! Start multinomial beta _ one value
      subroutine log_multinomial_beta_single(alpha,K,lik,max_iter,i)
        IMPLICIT NONE
        real*8, intent(in) :: alpha
        real*8,dimension(max_iter),intent(inout) :: lik
        integer*8,intent(in) :: K,max_iter,i

           lik(i) = lik(i) - (K * log_gamma(alpha)) - log_gamma(K * alpha)

      end
! End multinomial beta _ one value 

! Conditional Distribution

      subroutine conditional_distribution(matrix,NZW, NZM, NZ, beta, alpha,ntopics,M,N,p_z,j,ll)
      IMPLICIT NONE
        integer*8,intent(in) :: j,ll,M, N,ntopics
        integer*8 ii,nn,aa
        real*8, intent(inout) :: p_z(ntopics)
        real*8, intent(in) :: beta
        real*8, intent(in) :: alpha
        integer*8,dimension(N,M),intent(in) :: matrix
        integer*8,dimension(ntopics),intent(in) :: NZ
        integer*8,dimension(ntopics,M),intent(in) :: NZM
        integer*8,dimension(N,ntopics),intent(in) :: NZW
        integer*8 vsize
        
        vsize = 0
        do nn = 1,N
          if (matrix(nn,j) == 0) then
            cycle
          endif
          vsize = vsize + 1
        enddo
        do ii = 1,ntopics
          p_z(ii) = ((NZM(ii,j) + alpha) * (NZW(ll,ii) + beta)) / (NZ(ii) + vsize * beta)
        enddo
        ! this abs() shouldn't have to be here...
        p_z = abs(p_z)
        p_z = p_z / sum(p_z)
        
        ! We need to do something to remove 'bunching' to one topic
        ! since an error happens at genmul is any p_z is > .99
        ! for now, just iterate through, if a p_z is > .99 then
        ! subtract .01 from p_z(ii) and add .01 to p_z((ii-1))

        do ii = 1,ntopics
          
          if (0.99E+00 < p_z(ii)) then
            p_z(ii) = p_z(ii) - .1
            do aa = 1,ntopics
              if (ii == aa) then
                cycle
              endif
              p_z(aa) = p_z(aa) + (.1 / (ntopics-1))
            enddo
          elseif (.001 > p_z(ii)) then
            p_z(ii) = p_z(ii) + .1
            do aa = 1,ntopics
              if (ii == aa) then
                cycle
              endif
              p_z(aa) = p_z(aa) - (.1 / (ntopics-1))
            enddo
          endif
        enddo
      end

! End Conditional Distribution
END MODULE gibbs_sampler

