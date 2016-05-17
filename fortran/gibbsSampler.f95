! gibbsSampler.f95
MODULE gibbs_sampler

CONTAINS
      subroutine gibbsSampler(matrix,NZW,NZM,NZ,NM,ntopics, &
     max_iter,M,N,p_z,topics,alpha, beta,lik,top_size)
!!!!!!!
! Description: Main algorithm to run the Gibbs sampler for Latent 
!  Dirichlet Analysis.
!!!!!!!
! Everything must be row contiguous in fortran
! you can check this in python with <array>.flags.f_contiguous
        IMPLICIT NONE

        integer*8, dimension(n,m) :: matrix
        integer*8 m
        integer*8 n
        integer*8 top_size
        integer*8, dimension(ntopics,n) :: nzw
        integer*8, dimension(ntopics,m) :: nzm
        integer*8, dimension(ntopics) :: nz
        integer*8, dimension(m) :: nm
        integer*8, intent(in) :: ntopics
        integer*8, intent(in) :: max_iter
        integer*8, dimension(top_size) :: topics
        integer*8 :: i,j,ll,nn
        integer*8 :: Z, gg, kk
        integer*4, dimension(ntopics) :: genZ
        real*8, dimension(ntopics) :: p_z
        real*8,intent(in) :: alpha
        real*8,intent(in) :: beta
        real*8, dimension(max_iter), intent(inout) :: lik
        EXTERNAL genmul

      !kk is the iterator for topics
      do i=1,max_iter
       kk = 1
        do j=1,M
          do ll=1,N

            if (matrix(ll,j) == 0) then
              cycle
            endif

            do nn=1,matrix(ll,j)

              Z = topics(kk)
              NZM(Z,j) = NZM(Z,j) - 1
              NM(j) = NM(j) - 1
              NZW(Z,ll) = NZW(Z,ll) - 1
              NZ(Z) = NZ(Z) - 1

              call  conditional_distribution(matrix,NZW, NZM, NZ,NM, beta, alpha,ntopics,M,N,p_z,j,ll)

              genZ = 0
    
              call genmul(1,abs(p_z),ntopics,genZ)
              
              do gg = 1,ntopics
                if (genZ(gg) == 1) then
                  Z = gg
                  exit
                endif
              enddo

              topics(kk) = Z
              kk = kk + 1
              NZM(Z,j) = NZM(Z,j) + 1
              NM(j) = NM(j) + 1
              NZW(Z,ll) = NZW(Z,ll) + 1
              NZ(Z) = NZ(Z) + 1
            enddo  
          enddo
        enddo
        call loglikelihood(matrix,NZW, NZM, alpha, beta, ntopics,N,M,lik,max_iter,i,j)
        
        if (mod(i,10) == 0) then        
          write(*,*) 'Iteration:'
          write(*,*) i
          write(*,*) 'Likelihood:'
          write(*,*) lik(i)
        endif
      enddo
      end

! End gibbsSampler



! Begin Likelihood
      subroutine loglikelihood(matrix,NZW, NZM, alpha, beta, ntopics,N,M,lik,max_iter,i,j)
        IMPLICIT NONE
        integer*8,dimension(n,m),intent(in) :: matrix
        integer*8 vsize
        integer*8,intent(in) :: N, M, ntopics,max_iter,i,j
        integer*8 :: z,mm
        integer*8,dimension(ntopics,M),intent(in) :: NZM
        integer*8,dimension(ntopics,N),intent(in) :: NZW
        real*8, intent(in)  :: alpha, beta
        real*8,dimension(max_iter), intent(inout) :: lik
        
        vsize = count(matrix(:,j) /= 0)

        do z = 1,ntopics
              call log_multinomial_beta(NZW(:,Z) + beta,lik,ntopics,max_iter,i)
              ! Because below only takes in a single value we write seperate function
              call log_multinomial_beta_single(beta,vsize,lik,max_iter,i)

        enddo

        do mm = 1,M
              call log_multinomial_beta(NZM(:,mm) + alpha,lik,ntopics,max_iter,i)
              call log_multinomial_beta_single(alpha,ntopics,lik,max_iter,i)
        enddo

      end  
! End loglikelihood      

! Begin log_multinomial_beta
      subroutine log_multinomial_beta(alpha,lik,ntopics,max_iter,i)
        IMPLICIT NONE
        real*8,dimension(ntopics), intent(in) :: alpha
        real*8,dimension(max_iter),intent(inout) :: lik
        integer*8,intent(in) :: ntopics,max_iter,i

           lik(i) = lik(i) + sum(log_gamma(alpha)) - log_gamma(sum(alpha))
      end
     
! End multinomial beta

! Begin multinomial beta w/ one value
      subroutine log_multinomial_beta_single(alpha,K,lik,max_iter,i)
        IMPLICIT NONE
        real*8, intent(in) :: alpha
        real*8,dimension(max_iter),intent(inout) :: lik
        integer*8,intent(in) :: K,max_iter,i

           lik(i) = lik(i) - (K * log_gamma(alpha)) - log_gamma(K * alpha)

      end
! End multinomial beta w/ one value 


! Begin Conditional Distribution

      subroutine conditional_distribution(matrix,NZW, NZM, NZ,NM, beta, alpha,ntopics,M,N,p_z,j,ll)
      IMPLICIT NONE
        integer*8,intent(in) :: j,ll,M, N,ntopics
        integer*8 :: ii, aa
        real*8, intent(inout) :: p_z(ntopics)
        real*8, intent(in) :: beta
        real*8, intent(in) :: alpha
        real*8,dimension(ntopics) :: left
        real*8,dimension(ntopics) :: right
        real*8 :: tolup, toldown
        real*8 shift
        integer*8, dimension(M) :: NM
        integer*8,dimension(N,M),intent(in) :: matrix
        integer*8,dimension(ntopics),intent(in) :: NZ
        integer*8,dimension(ntopics,M),intent(in) :: NZM
        integer*8,dimension(N,ntopics),intent(in) :: NZW
        integer*8 vsize
        
        tolup = .999999999999999
        toldown = .0000000000000001
        shift = .000000000000001
        

        vsize = count(matrix(:,j) /= 0)

        left = (NZW(ll,:) + beta) / (NZ + beta * vsize)

        right = (NZM(:,j) + alpha) / (NM + alpha * ntopics)

        p_z = left * right
        
        p_z = p_z / sum(p_z)
        ! We need to do something to remove 'bunching' to one topic
        ! since an error happens at genmul is any p_z is > .99
        ! for now, just iterate through, if a p_z is > .99 then
        ! subtract .01 from p_z(ii) and add .01 to p_z((ii-1))


        do ii = 1,ntopics
          
          if (tolup < p_z(ii)) then
            p_z(ii) = p_z(ii) - shift
            do aa = 1,ntopics
              if (ii == aa) then
                cycle
              endif
              write(*,*) 'WARNING: Tol adjusted for > tol '
              p_z(aa) = p_z(aa) + (shift / (ntopics-1))
            enddo
          elseif (toldown > p_z(ii)) then
            p_z(ii) = p_z(ii) + shift
            do aa = 1,ntopics
              if (ii == aa) then
                cycle
              endif
              write(*,*) 'WARNING: Tol adjusted for < tol'
              p_z(aa) = p_z(aa) - (shift / (ntopics-1))
            enddo
          endif
        enddo
        p_z = p_z / sum(p_z)

      end
! End Conditional Distribution
END MODULE gibbs_sampler

