!#####################################################################
!# wrapper functions for norm2 shared object
!#####################################################################
subroutine norm_em(n, r, p, x, y, mvcode, &
        prior_type_int, prior_df, prior_sscp, &
        max_iter, criterion, estimate_worst, &
        startval_present, &
        iter, converged, reldiff, loglik, logpost, &
        beta, sigma, yimp, &
        npatt, mis, n_in_patt, nobs, which_patt, ybar, ysdv, &
        rate_beta, rate_sigma, em_worst_ok, worst_frac, &
        nparam, worst_linear_coef, &
        status, msg_len_max, msg_codes, msg_len_actual)
   !DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:"norm_em_" :: norm_em
   !#############################################################
   ! This is a wrapper function for run_norm_engine_em
   ! Values for prior_type_int:
   !    1 : uniform
   !    2 : jeffreys
   !    3 : ridge (prior_df is relevant)
   !    4 : invwish (prior.df and prior.sscp are relevant)
   ! If startval_present is .TRUE., input values of beta and sigma
   ! are interpreted as starting values; otherwise, defaults are
   ! used.
   ! status = 0 means everything ran OK, or EM was aborted for
   !    some reason (see msg)
   ! Other value of status indicates a fatal error.
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use norm_engine
   use random_generator
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n  ! # rows of y
   integer(kind=our_int), intent(in) :: r  ! # columns of y
   integer(kind=our_int), intent(in) :: p  ! # columns of x
   real(kind=our_dble), intent(in) :: x(n,p) ! predictor variables
   real(kind=our_dble), intent(in) :: y(n,r) ! response variables
   real(kind=our_dble), intent(in) :: mvcode ! missing value code
   integer(kind=our_int), intent(in) :: prior_type_int
   real(kind=our_dble), intent(in) :: prior_df
   real(kind=our_dble), intent(in) :: prior_sscp(r,r)
   integer(kind=our_int), intent(in) :: max_iter
   real(kind=our_dble), intent(in) :: criterion
   logical, intent(in) :: estimate_worst
   logical, intent(in) :: startval_present
   ! declare output arguments
   integer(kind=our_int), intent(out) :: iter
   logical, intent(out) :: converged
   real(kind=our_dble), intent(out) :: reldiff
   real(kind=our_dble), intent(out) :: loglik(max_iter)
   !   Note: the elements of loglik beyond "iter" are irrelevant.
   real(kind=our_dble), intent(out) :: logpost(max_iter)
   !   Note: the elements of logpost beyond "iter" are irrelevant.
   real(kind=our_dble), intent(inout) :: beta(p,r)
   real(kind=our_dble), intent(inout) :: sigma(r,r)
   real(kind=our_dble), intent(out) :: yimp(n,r)
   integer(kind=our_int), intent(out) :: npatt
   logical, intent(out) :: mis(n,r)
   !   Note: the rows of mis beyond "npatt" are irrelevant.
   integer(kind=our_int), intent(out) :: n_in_patt(n)
   !   Note: the elements of n_in_patt beyond "npatt" are irrelevant.
   integer(kind=our_int), intent(out) :: nobs(r)
   integer(kind=our_int), intent(out) :: which_patt(n)
   real(kind=our_dble), intent(out) :: ybar(r)
   real(kind=our_dble), intent(out) :: ysdv(r)
   real(kind=our_dble), intent(out) :: rate_beta(p,r)
   real(kind=our_dble), intent(out) :: rate_sigma(r,r)
   logical, intent(out) :: em_worst_ok
   real(kind=our_dble), intent(out) :: worst_frac
   integer(kind=our_int), intent(in) :: nparam
   real(kind=our_dble), intent(out) :: worst_linear_coef(nparam)
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 8 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare local pointers to pass as arguments
   real(kind=our_dble), pointer :: beta_start_local(:,:), &
        sigma_start_local(:,:), loglik_local(:), logpost_local(:), &
        beta_local(:,:), sigma_local(:,:), yimp_local(:,:)
   logical, pointer :: mis_local(:,:)
   integer(kind=our_int), pointer :: n_in_patt_local(:), &
        nobs_local(:), which_patt_local(:)
   real(kind=our_dble), pointer :: &
        ybar_local(:), ysdv_local(:), &
        rate_beta_local(:,:), rate_sigma_local(:,:), &
        worst_linear_coef_local(:)
   ! other locals
   character(len=20) :: prior_type
   integer(our_int) :: ijunk
   character(len=*), parameter :: platform = "S+"
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( dyn_dealloc(   beta_start_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  sigma_start_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( loglik_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( logpost_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(   beta_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  sigma_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  yimp_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  mis_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  n_in_patt_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  nobs_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  which_patt_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  ybar_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  ysdv_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  rate_beta_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc(  rate_sigma_local, err ) == RETURN_FAIL ) goto 800
   ! allocate local pointers that need it
   if( dyn_alloc( beta_start_local, p, r, err ) == RETURN_FAIL ) goto 800
   if( dyn_alloc( sigma_start_local, r, r, err ) == RETURN_FAIL ) goto 800
   if( startval_present) then
      beta_start_local(:,:) = beta(:,:)
      sigma_start_local(:,:) = sigma(:,:)
   end if
   ! set prior_type
   if( prior_type_int == 1 ) then
      prior_type = "uniform"
   else if( prior_type_int == 2 ) then
      prior_type = "jeffreys"
   else if( prior_type_int == 3 ) then
      prior_type = "ridge"
   else if( prior_type_int == 4 ) then
      prior_type = "invwish"
   else
      prior_type = "other"
   end if
   if( run_norm_engine_em(x = x, y = y, mvcode = mvcode, &
        startval_present = startval_present, &
        beta_start = beta_start_local, &
        sigma_start = sigma_start_local, &
        prior_type = trim(prior_type), &
        prior_df = prior_df, &
        prior_sscp = prior_sscp, &
        max_iter = max_iter, &
        criterion = criterion, &
        estimate_worst = estimate_worst, &
        err = err, &
        iter = iter, &
        converged = converged, &
        reldiff = reldiff, &
        loglik = loglik_local, &
        logpost = logpost_local, &
        beta = beta_local, &
        sigma = sigma_local, &
        yimp = yimp_local, &
        mis = mis_local, &
        n_in_patt = n_in_patt_local, &
        nobs = nobs_local, &
        which_patt = which_patt_local, &
        ybar = ybar_local, &
        ysdv = ysdv_local, &
        rate_beta = rate_beta_local, &
        rate_sigma = rate_sigma_local, &
        em_worst_ok = em_worst_ok, &
        worst_frac = worst_frac, &
        worst_linear_coef = worst_linear_coef_local &
        ) == RETURN_FAIL ) goto 800
   if( associated( loglik_local ) ) then
      loglik(:) = 0.D0
      loglik(1:iter) = loglik_local(1:iter)
   end if
   if( associated( logpost_local ) )  then
      logpost(:) = 0.D0
      logpost(1:iter) = logpost_local(1:iter)
   end if
   if( associated( beta_local ) ) &
        beta(:,:) = beta_local(:,:)
   if( associated( sigma_local ) ) &
        sigma(:,:) = sigma_local(:,:)
   if( associated( yimp_local ) ) &
        yimp(:,:) = yimp_local(:,:)
   if( associated( mis_local ) ) then
      npatt = size( mis_local, 1 )
      mis(:,:) = .false.
      mis( 1:npatt, 1:r ) = mis_local(:,:)
   end if
   if( associated( n_in_patt_local ) ) then
      n_in_patt(:) = 0
      n_in_patt( 1:npatt ) = n_in_patt_local(:)
   end if
   if( associated( nobs_local ) ) &
      nobs(:) = nobs_local(:)
   if( associated( which_patt_local ) ) &
        which_patt(:) = which_patt_local(:)
   if( associated( ybar_local ) ) &
        ybar(:) = ybar_local(:)
   if( associated( ysdv_local ) ) &
        ysdv(:) = ysdv_local(:)
   if( associated( rate_beta_local ) ) &
        rate_beta(:,:) = rate_beta_local(:,:)
   if( associated( rate_sigma_local ) ) &
        rate_sigma(:,:) = rate_sigma_local(:,:)
   if( associated( worst_linear_coef_local ) ) &
        worst_linear_coef(:) = worst_linear_coef_local(:)
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = dyn_dealloc( beta_start_local, err )
   ijunk = dyn_dealloc( sigma_start_local, err )
   ijunk = dyn_dealloc( loglik_local, err )
   ijunk = dyn_dealloc( logpost_local, err )
   ijunk = dyn_dealloc( beta_local, err )
   ijunk = dyn_dealloc( sigma_local, err )
   ijunk = dyn_dealloc( yimp_local, err )
   ijunk = dyn_dealloc( mis_local, err )
   ijunk = dyn_dealloc( n_in_patt_local, err )
   ijunk = dyn_dealloc( nobs_local, err )
   ijunk = dyn_dealloc( which_patt_local, err )
   ijunk = dyn_dealloc( ybar_local, err )
   ijunk = dyn_dealloc( ysdv_local, err )
   ijunk = dyn_dealloc( rate_beta_local, err )
   ijunk = dyn_dealloc( rate_sigma_local, err )
   ijunk = dyn_dealloc( worst_linear_coef_local, err )
end subroutine norm_em
!#####################################################################
subroutine norm_mcmc(n, r, p, x, y, mvcode, &
        prior_type_int, prior_df, prior_sscp, &
        iter, multicycle, seeds, &
        impute_every, nimps, &
        save_all_series, save_worst_series, worst_linear_coef, &
        series_length, &
        beta, sigma, yimp,  iter_actual, beta_series, sigma_series, &
        loglik, logpost, worst_series, &
        npatt, mis, n_in_patt, nobs, which_patt, ybar, ysdv, &
        imp_list, nimps_actual, &
        status, msg_len_max, msg_codes, msg_len_actual)
   !DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:"norm_mcmc_" :: norm_mcmc
   !#############################################################
   ! This is a wrapper function for run_norm_engine_mcmc
   ! Values for prior_type_int:
   !    1 : uniform
   !    2 : jeffreys
   !    3 : ridge (prior_df is relevant)
   !    4 : invwish (prior.df and prior.sscp are relevant)
   ! If multiple imputations are desired, then argument nimps 
   ! (which is also the last dimension of imp_list) should
   ! be set to floor( iter / impute_every ).
   ! If no multiple imputations are desired, then impute_every
   ! and nimps should be set to zero.
   ! If save_all_series is .true., then series_length should be set 
   !  equal to iter.  Otherwise, set series_length to zero.   
   ! status = 0 means everything ran OK, or MCMC was aborted for
   !    some reason (see msg)
   ! iter_actual is the actual number of iterations
   !   performed (will be equal to iter, unless aborted)
   ! nimps_actual is the actual number of imputations
   !   created (will be equal to nimps, unless aborted)
   ! Other value of status indicates a fatal error
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use norm_engine
   use random_generator
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n  ! # rows of y
   integer(kind=our_int), intent(in) :: r  ! # columns of y
   integer(kind=our_int), intent(in) :: p  ! # columns of x
   real(kind=our_dble), intent(in) :: x(n,p) ! predictor variables
   real(kind=our_dble), intent(in) :: y(n,r) ! response variables
   real(kind=our_dble), intent(in) :: mvcode ! missing value code
   integer(kind=our_int), intent(in) :: prior_type_int
   real(kind=our_dble), intent(in) :: prior_df
   real(kind=our_dble), intent(in) :: prior_sscp(r,r)
   integer(kind=our_int), intent(in) :: iter
   integer(kind=our_int), intent(in) :: multicycle
   integer(kind=our_int), intent(in) :: seeds(2)
   integer(kind=our_int), intent(in) :: impute_every
   integer(kind=our_int), intent(in) :: nimps
   logical, intent(in) :: save_all_series, save_worst_series
   real(kind=our_dble), intent(in) :: &
        worst_linear_coef( p*r + r*(r+1)/2 )
   integer(kind=our_int), intent(in) :: series_length
   ! declare inout arguments
   real(kind=our_dble), intent(inout) :: beta(p,r)
   real(kind=our_dble), intent(inout) :: sigma(r,r)
   ! declare output arguments
   real(kind=our_dble), intent(out) :: yimp(n,r)
   integer(kind=our_int), intent(out) :: iter_actual
   !  Note: the elements of these series beyond
   !   "iter_actual" are irrelevant
   real(kind=our_dble), intent(out) :: beta_series(p,r,series_length)
   real(kind=our_dble), intent(out) :: sigma_series(r,r,series_length)
   real(kind=our_dble), intent(out) :: loglik(iter)
   real(kind=our_dble), intent(out) :: logpost(iter)
   real(kind=our_dble), intent(out) :: worst_series(iter)
   integer(kind=our_int), intent(out) :: npatt
   logical, intent(out) :: mis(n,r)
   !   Note: the rows of mis beyond "npatt" are irrelevant.
   integer(kind=our_int), intent(out) :: n_in_patt(n)
   !   Note: the elements of n_in_patt beyond "npatt" are irrelevant.
   integer(kind=our_int), intent(out) :: nobs(r)
   integer(kind=our_int), intent(out) :: which_patt(n)
   real(kind=our_dble), intent(out) :: ybar(r)
   real(kind=our_dble), intent(out) :: ysdv(r)
   real(kind=our_dble), intent(out) :: imp_list(n,r,nimps)
   !  Note: the elements of imp_list beyond "nimps_actual" are irrelevant
   integer(kind=our_int), intent(out) :: nimps_actual
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 8 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   character(len=20) :: prior_type
   real(kind=our_dble), pointer :: beta_start_local(:,:), &
        sigma_start_local(:,:), beta_local(:,:), &
        sigma_local(:,:), yimp_local(:,:), &
        beta_series_local(:,:,:), sigma_series_local(:,:,:), &
        worst_series_local(:), &
        loglik_local(:), logpost_local(:), &
        ybar_local(:), ysdv_local(:), &
        imp_list_local(:,:,:)
   logical, pointer :: mis_local(:,:)
   integer(kind=our_int), pointer :: n_in_patt_local(:), &
        nobs_local(:), which_patt_local(:)
   type(random_gendata) :: rand
   integer(our_int) :: iseed1, iseed2, ijunk
   character(len=*), parameter :: platform = "S+"
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( dyn_dealloc( beta_start_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( sigma_start_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( beta_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( sigma_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( yimp_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( beta_series_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( sigma_series_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( worst_series_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( loglik_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( logpost_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( ybar_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( ysdv_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( imp_list_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( n_in_patt_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( nobs_local, err ) == RETURN_FAIL ) goto 800
   if( dyn_dealloc( which_patt_local, err ) == RETURN_FAIL ) goto 800
   ! allocate local pointers that need it
   if( dyn_alloc( beta_start_local, p, r, err ) == RETURN_FAIL ) goto 800
   if( dyn_alloc( sigma_start_local, r, r, err ) == RETURN_FAIL ) goto 800
   beta_start_local(:,:) = beta(:,:)
   sigma_start_local(:,:) = sigma(:,:)
   ! set prior_type
   if( prior_type_int == 1 ) then
      prior_type = "uniform"
   else if( prior_type_int == 2 ) then
      prior_type = "jeffreys"
   else if( prior_type_int == 3 ) then
      prior_type = "ridge"
   else if( prior_type_int == 4 ) then
      prior_type = "invwish"
   else
      prior_type = "other"
   end if
   ! set random seeds
   iseed1 = seeds(1)
   iseed2 = seeds(2)
   if( ran_setall( rand, iseed1, iseed2, err) == RETURN_FAIL ) goto 800
   ! call fortran
   if( run_norm_engine_mcmc(x = x, y = y, &
        mvcode = mvcode, &
        rand = rand, &
        beta_start = beta_start_local, & 
        sigma_start = sigma_start_local, &
        prior_type = trim(prior_type), &
        prior_df = prior_df, &
        prior_sscp = prior_sscp, &
        max_iter = iter, &
        multicycle = multicycle, &
        impute_every = impute_every, &
        save_all_series = save_all_series, &
        save_worst_series = save_worst_series, &
        worst_linear_coef = worst_linear_coef, &
        err = err, &
        iter = iter_actual, &
        beta = beta_local, &
        sigma = sigma_local, &
        loglik = loglik_local, &
        logpost = logpost_local, &
        yimp = yimp_local, &
        mis = mis_local, &
        n_in_patt = n_in_patt_local, &
        nobs = nobs_local, &
        which_patt = which_patt_local, &
        ybar = ybar_local, &
        ysdv = ysdv_local, &
        beta_series = beta_series_local, & ! if save_all_series
        sigma_series = sigma_series_local, & ! if save_all_series
        worst_series = worst_series_local, &
        nimp = nimps_actual, &
        imp_list = imp_list_local ) == RETURN_FAIL ) goto 800
   if( associated( beta_local ) ) &
        beta(:,:) = beta_local(:,:)
   if( associated( sigma_local ) ) &
        sigma(:,:) = sigma_local(:,:)
   if( associated( loglik_local ) ) then
      loglik(:) = 0.D0
      loglik(1:iter_actual) = loglik_local(1:iter_actual)
   end if
   if( associated( logpost_local ) ) then
      logpost(:) = 0.D0
      logpost(1:iter_actual) = logpost_local(1:iter_actual)
   end if
   if( associated( yimp_local ) ) &
        yimp(:,:) = yimp_local(:,:)
   if( associated( mis_local ) ) then
      npatt = size( mis_local, 1 )
      mis(:,:) = .false.
      mis( 1:npatt, 1:r ) = mis_local(:,:)
   end if
   if( associated( n_in_patt_local ) ) then
      n_in_patt(:) = 0
      n_in_patt( 1:npatt ) = n_in_patt_local(:)
   end if
   if( associated( nobs_local ) ) &
        nobs(:) = nobs_local(:)
   if( associated( which_patt_local ) ) &
        which_patt(:) = which_patt_local(:)
   if( associated( ybar_local ) ) &
        ybar(:) = ybar_local(:)
   if( associated( ysdv_local ) ) &
        ysdv(:) = ysdv_local(:)
   if( save_all_series ) then
      beta_series(:,:,:) = 0.D0
      beta_series(:,:,1:iter_actual) = &
           beta_series_local(:,:,1:iter_actual)
      sigma_series(:,:,:) = 0.D0
      sigma_series(:,:,1:iter_actual) &
           = sigma_series_local(:,:,1:iter_actual)
   end if
   if( save_worst_series .and. associated( worst_series_local ) ) then
      worst_series(:) = 0.D0
      worst_series(1:iter_actual) = worst_series_local(1:iter_actual)
   else
      worst_series(:) = 0.D0
   end if
   if( impute_every /= 0 ) then
      imp_list(:,:,:) = 0.D0
      imp_list(:,:,1:nimps_actual) = imp_list_local(:,:,1:nimps_actual)
   end if
  ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = dyn_dealloc( beta_start_local, err)
   ijunk = dyn_dealloc( sigma_start_local, err)
   ijunk = dyn_dealloc( beta_local, err)
   ijunk = dyn_dealloc( sigma_local, err)
   ijunk = dyn_dealloc( yimp_local, err)
   ijunk = dyn_dealloc( beta_series_local, err)
   ijunk = dyn_dealloc( sigma_series_local, err)
   ijunk = dyn_dealloc( worst_series_local, err)
   ijunk = dyn_dealloc( loglik_local, err)
   ijunk = dyn_dealloc( logpost_local, err)
   ijunk = dyn_dealloc( ybar_local, err)
   ijunk = dyn_dealloc( ysdv_local, err)
   ijunk = dyn_dealloc( imp_list_local, err)
   ijunk = dyn_dealloc( n_in_patt_local, err)
   ijunk = dyn_dealloc( nobs_local, err)
   ijunk = dyn_dealloc( which_patt_local, err)
end subroutine norm_mcmc
!#####################################################################
subroutine norm_imp_rand(n, r, p, x, y, mvcode, seeds, &
        beta, sigma, yimp, status, &
        msg_len_max, msg_codes, msg_len_actual)
   !DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:"norm_imp_rand_" :: norm_imp_rand
   !#############################################################
   ! This is a wrapper function for run_norm_engine_impute_random
   ! status = 0 means everything ran OK
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use norm_engine
   use random_generator
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n  ! # rows of y
   integer(kind=our_int), intent(in) :: r  ! # columns of y
   integer(kind=our_int), intent(in) :: p  ! # columns of x
   real(kind=our_dble), intent(in) :: x(n,p) ! predictor variables
   real(kind=our_dble), intent(in) :: y(n,r) ! response variables
   real(kind=our_dble), intent(in) :: mvcode ! missing value code
   integer(kind=our_int), intent(in) :: seeds(2)
   real(kind=our_dble), intent(in) :: beta(p,r)
   real(kind=our_dble), intent(in) :: sigma(r,r)
   ! declare output arguments
   real(kind=our_dble), intent(out) :: yimp(n,r)
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 8 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   real(kind=our_dble), pointer :: yimp_local(:,:)
   type(random_gendata) :: rand
   integer(our_int) :: iseed1, iseed2, ijunk
   character(len=*), parameter :: platform = "S+"
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( dyn_dealloc( yimp_local, err ) == RETURN_FAIL ) goto 800
   iseed1 = seeds(1)
   iseed2 = seeds(2)
   if( ran_setall( rand, iseed1, iseed2, err) == RETURN_FAIL ) goto 800
   if( run_norm_engine_impute_random(x = x, y = y, &
           mvcode = mvcode, &
           beta = beta, &
           sigma = sigma, &
           yimp = yimp_local, &
           rand = rand, &
           err = err) == RETURN_FAIL ) goto 800
   if( associated( yimp_local ) ) then
      yimp(:,:) = yimp_local(:,:)
   else
      goto 800
   end if
  ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = dyn_dealloc( yimp_local, err)
end subroutine norm_imp_rand
!#####################################################################
subroutine norm_imp_mean(n, r, p, x, y, mvcode, &
        beta, sigma, yimp, status, &
        msg_len_max, msg_codes, msg_len_actual)
   !DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:"norm_imp_mean_" :: norm_imp_mean
   !#############################################################
   ! This is a wrapper function for run_norm_engine_impute_mean
   ! status = 0 means everything ran OK
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use norm_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n  ! # rows of y
   integer(kind=our_int), intent(in) :: r  ! # columns of y
   integer(kind=our_int), intent(in) :: p  ! # columns of x
   real(kind=our_dble), intent(in) :: x(n,p) ! predictor variables
   real(kind=our_dble), intent(in) :: y(n,r) ! response variables
   real(kind=our_dble), intent(in) :: mvcode ! missing value code
   real(kind=our_dble), intent(in) :: beta(p,r)
   real(kind=our_dble), intent(in) :: sigma(r,r)
   ! declare output arguments
   real(kind=our_dble), intent(out) :: yimp(n,r)
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 8 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   real(kind=our_dble), pointer :: yimp_local(:,:)
   integer(our_int) :: ijunk
   character(len=*), parameter :: platform = "S+"
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   if( dyn_dealloc( yimp_local, err ) == RETURN_FAIL ) goto 800
   if( run_norm_engine_impute_mean(x = x, y = y, &
           mvcode = mvcode, &
           beta = beta, &
           sigma = sigma, &
           yimp = yimp_local, &
           err = err) == RETURN_FAIL ) goto 800
   if( associated( yimp_local ) ) then
      yimp(:,:) = yimp_local(:,:)
   else
      goto 800
   end if
  ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
   ijunk = dyn_dealloc( yimp_local, err)
end subroutine norm_imp_mean
!#####################################################################
subroutine norm_logpost(n, r, p, x, y, mvcode, &
        beta, sigma, prior_type_int, prior_df, prior_sscp, &
        logpost, &
        status, msg_len_max, msg_codes, msg_len_actual)
   !DEC$ ATTRIBUTES DLLEXPORT,C,REFERENCE,ALIAS:"norm_logpost_" :: norm_logpost
   !#############################################################
   ! This is a wrapper function for run_norm_engine_logpost
   ! Values for prior_type_int:
   !    1 : uniform
   !    2 : jeffreys
   !    3 : ridge (prior_df is relevant)
   !    4 : invwish (prior.df and prior.sscp are relevant)
   ! status = 0 means everything ran OK
   !#############################################################
   use error_handler
   use program_constants
   use dynalloc
   use norm_engine
   implicit none
   ! declare input arguments
   integer(kind=our_int), intent(in) :: n  ! # rows of y
   integer(kind=our_int), intent(in) :: r  ! # columns of y
   integer(kind=our_int), intent(in) :: p  ! # columns of x
   real(kind=our_dble), intent(in) :: x(n,p) ! predictor variables
   real(kind=our_dble), intent(in) :: y(n,r) ! response variables
   real(kind=our_dble), intent(in) :: mvcode ! missing value code
   real(kind=our_dble), intent(in) :: beta(p,r)
   real(kind=our_dble), intent(in) :: sigma(r,r)
   integer(kind=our_int), intent(in) :: prior_type_int
   real(kind=our_dble), intent(in) :: prior_df
   real(kind=our_dble), intent(in) :: prior_sscp(r,r)
   ! declare output arguments
   real(kind=our_dble), intent(out) :: logpost
   integer(kind=our_int), intent(out) :: status
   integer(kind=our_int), intent(in) :: msg_len_max
   integer(kind=our_int), intent(out) :: msg_codes( msg_len_max, 8 )
   integer(kind=our_int), intent(out) :: msg_len_actual
   ! declare locals
   character(len=20) :: prior_type
   character(len=*), parameter :: platform = "S+"
   type(error_type) :: err
   ! begin
   status = 1
   call err_reset(err)
   ! set prior_type
   if( prior_type_int == 1 ) then
      prior_type = "uniform"
   else if( prior_type_int == 2 ) then
      prior_type = "jeffreys"
   else if( prior_type_int == 3 ) then
      prior_type = "ridge"
   else if( prior_type_int == 4 ) then
      prior_type = "invwish"
   else
      prior_type = "other"
   end if
   if( run_norm_engine_logpost(x = x, y = y, &
        mvcode = mvcode, &
        beta = beta, &
        sigma = sigma, &
        logpost = logpost, &
        err = err, &
        prior_type = trim(prior_type), &
        prior_df = prior_df, &
        prior_sscp = prior_sscp) == RETURN_FAIL ) goto 800
   ! normal exit
   status = 0
800 continue
   ! report message if present
   msg_codes(:,:) = 0
   msg_len_actual = 0
   if( err_msg_present(err) ) call err_get_codes(err, &
        msg_codes, msg_len_actual)
   ! cleanup
   call err_reset(err)
end subroutine norm_logpost
!#####################################################################
