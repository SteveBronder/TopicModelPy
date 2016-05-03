! gibbsSampler.f95
MODULE gibbs_sampler

CONTAINS
      subroutine gibbsSampler(matrix,NZW,NZM,NZ,NM,ntopics, &
     max_iter,M,N,p_z,topics,topics2)
        IMPLICIT NONE
! Everything is transpose in fortran (??) so rows and columns get switched
        integer :: M, N,i,j,ll,nn,ntapp
        integer max_iter
        integer Z
        integer matrix(M,N)
        integer ntopics
        real*4 p_z(ntopics)
        integer topics(m,n)
        integer topics2(m,n)
        integer NZ(ntopics)
        integer NM(M)
        integer NZM(ntopics,M)
        integer NZW(N,ntopics)
        real*4 :: alpha, beta
        EXTERNAL genmul
      ntapp = ntopics-1
      do i=1,max_iter
        if (i > 1) then
          topics = topics2
        endif
        do j=1,M
! FIXME: ll is like w in the enumerate
!   (i.e. if first nonzero is at place 27, rep 27 (ll) 4 times (ll)
          do ll=1,N
            if (matrix(j,ll) == 0) then
              cycle
            endif
            ! note: j = M, ll = N, nn = w
            do nn=1,matrix(j,ll)
              Z = topics(m,n)
              ! Note: Due to memory access in fortran columns have to be rows
              !  This is why these are reversed and the transpose is
              !  brought in
              NZM(Z,j) = NZM(Z,j) - 1
              NM(j) = NM(j) - 1
              NZW(n,Z) = NZW(n,Z) - 1
              NZ(Z) = NZ(Z) - 1

              call  conditional_distribution(matrix,NZW, NZM, NZ, NM, beta, alpha,ntopics,M,N,p_z,j,ll)
              call genmul(1,abs(p_z(1:ntapp)),ntopics,Z)
              topics2(m,n) = Z
              NZM(Z,j) = NZM(Z,j) + 1
              NM(j) = NM(j) + 1
              NZW(n,Z) = NZW(n,Z) + 1
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

      subroutine conditional_distribution(matrix,NZW, NZM, NZ, NM, beta, alpha,ntopics,M,N,p_z,j,ll)
      IMPLICIT NONE
        integer :: j,ll, M, N,i
        integer ntopics
        real*4 p_z(ntopics)
        real*4 beta
        real*4 alpha
        real*4 left(ntopics,(N-1))
        real*4 right(N-1)
        integer matrix(M,N)
        integer matrix_cut(M,(N-1))
        integer NZ(ntopics)
        integer NM(M)
        integer NZM(ntopics,M)
        integer NZW(N,ntopics)
        integer NZW_cut((N-1),ntopics) 

        ! ll is the row to exclude
        do i=1,N
          if (i == ll) then
            cycle
          endif
          NZW_cut(i,:) = NZW(i,:)
        enddo

        do i=1,N
          if (i == ll) then
            cycle
          endif
          matrix_cut(:,i) = matrix(:,i)
        enddo

        ! Pretty sure this should not be N, but V
        ! V:= number of unique words in that document
        left = (NZW_cut + beta) / sum(NZ + beta)
        right = (matrix_cut(j,:) + alpha) / (NM + alpha)
        p_z = matmul(left , right(1:(N-1)))
        p_z = abs(p_z)
        p_z = p_z / sum(p_z)
      end
! End Conditional Distribution
END MODULE gibbs_sampler

