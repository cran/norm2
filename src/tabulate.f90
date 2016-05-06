!#######################################################################
module tabulate
   ! Module for tabulating double precision, integer or
   ! character string variables using binary trees
   use program_constants
   use dynalloc
   use error_handler
   implicit none
   private
   public :: table_type  ! public type for storing the table
   public :: nullify_table, tabulate_variable, get_table_length, &
        get_table_type, get_table_values, get_table_frequencies
   interface tabulate_variable
      module procedure tabulate_real_variable
      module procedure tabulate_integer_variable
      module procedure tabulate_string_variable
   end interface
   interface get_table_values
      module procedure get_real_table_values
      module procedure get_integer_table_values
      module procedure get_string_table_values
   end interface
   !##################################################################
   integer(our_int), parameter :: max_string_length=132
   character(len=*), parameter :: modname = "tabulate"
   !##################################################################
   type dnode
      ! private type
      sequence
      real(kind=our_dble) :: value = 0.D0
      integer(kind=our_int) :: frequency = 0
      type(dnode), pointer :: left=>null()
      type(dnode), pointer :: right=>null()
   end type dnode
   !##################################################################
   type inode
      ! private type
      sequence
      integer(kind=our_int) :: value = 0
      integer(kind=our_int) :: frequency = 0
      type(inode), pointer :: left=>null()
      type(inode), pointer :: right=>null()
   end type inode
   !##################################################################
   type snode
      ! private type
      sequence
      character(len=max_string_length) :: value = ""
      integer(kind=our_int) :: frequency = 0
      type(snode), pointer :: left=>null()
      type(snode), pointer :: right=>null()
   end type snode
   !##################################################################
   type :: table_type
      ! public type whose contents are private
      private
      sequence
      logical :: is_null = .true.
      character(len=20) :: type = ""
      integer(kind=our_int) :: length=0
      integer(kind=our_int), pointer :: freq(:)=>null()
      real(kind=our_dble), pointer :: rVal(:)=>null()
      integer(kind=our_int), pointer :: iVal(:)=>null()
      character(len=max_string_length), pointer :: sVal(:)=>null()
      integer(kind=our_int) :: write_posn=0
   end type table_type
   !##################################################################
contains
   !################################################################
   integer(kind=our_int) function nullify_table(table, err) &
        result(answer)
      ! Nullifies a table_type object
      implicit none
      ! declare args
      type(table_type), intent(inout) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "nullify_table"
      ! go
      answer = RETURN_FAIL
      table%is_null = .false.
      table%type = ""
      table%length = 0
      table%write_posn = 0
      if( dyn_dealloc(table%freq, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(table%rVal, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(table%iVal, err) == RETURN_FAIL ) goto 800
      if( dyn_dealloc(table%sVal, err) == RETURN_FAIL ) goto 800
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function nullify_table
   !################################################################
   integer(kind=our_int) function tabulate_real_variable(vec, &
        table, err) result(answer)
      ! tabulates a dble vector.
      implicit none
      ! args
      real(kind=our_dble), intent(in) :: vec(:)
      type(table_type), intent(out) :: table
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i, ijunk
      real(kind=our_dble) :: number
      character(len=12) :: sInt
      type(dnode), pointer :: tree
      character(len=*), parameter :: subname = "tabulate_real_variable"
      ! begin
      answer = RETURN_FAIL
      nullify(tree)
      if( nullify_table(table,err) == RETURN_FAIL ) goto 800
      do i = 1, size(vec)
         number = vec(i)
         if( tabulate_real_number( tree, number, err) &
              == RETURN_FAIL ) goto 700
      end do
      if( write_real_tree_to_table(tree, table, err) == RETURN_FAIL ) &
           goto 800
      if( kill_real_tree(tree, err) == RETURN_FAIL ) goto 800
      ! normal exit
      table%is_null = .false.
      table%type = "double precision"
      answer = RETURN_SUCCESS
      return
      ! error traps
700   write(sInt,"(I12)") i
      sInt = adjustl(sInt)
      call err_handle(err, 1, &
           comment = "Failed to tabulate data" )
      call err_handle(err, 3, iobs = i )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
   end function tabulate_real_variable
   !################################################################
   integer(kind=our_int) function tabulate_integer_variable(vec, &
        table, err) result(answer)
      ! tabulates an integer vector.
      implicit none
      ! args
      integer(kind=our_int), intent(in) :: vec(:)
      type(table_type), intent(out) :: table
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i, ijunk
      integer(kind=our_int) :: number
      character(len=12) :: sInt
      type(inode), pointer :: tree
      character(len=*), parameter :: subname = "tabulate_integer_variable"
      ! begin
      answer = RETURN_FAIL
      nullify(tree)
      do i = 1, size(vec)
         number = vec(i)
         if( tabulate_integer_number( tree, number, err) &
              == RETURN_FAIL ) goto 700
      end do
      if( write_integer_tree_to_table(tree, table, err) == RETURN_FAIL ) &
           goto 800
      if( kill_integer_tree(tree, err) == RETURN_FAIL ) goto 800
      ! normal exit
      table%is_null = .false.
      table%type = "integer"
      answer = RETURN_SUCCESS
      return
      ! error traps
700   write(sInt,"(I12)") i
      sInt = adjustl(sInt)
      call err_handle(err, 1, &
           comment = "Failed to tabulate data" )
      call err_handle(err, 3, iobs = i )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
   end function tabulate_integer_variable
   !################################################################
   integer(kind=our_int) function tabulate_string_variable(vec, &
        table, err) result(answer)
      ! tabulates a character vector.
      implicit none
      ! args
      character(len=*), intent(in) :: vec(:)
      type(table_type), intent(out) :: table
      type(error_type), intent(inout) :: err
      ! locals
      integer(kind=our_int) :: i, ijunk
      character(len=max_string_length) :: string
      character(len=12) :: sInt
      type(snode), pointer :: tree
      character(len=*), parameter :: subname = "tabulate_string_variable"
      ! begin
      answer = RETURN_FAIL
      if( size( vec ) == 0 ) goto 600
      if( len( vec(1) ) > max_string_length ) goto 650
      nullify(tree)
      do i = 1, size(vec)
         string = vec(i)
         if( tabulate_string( tree, string, err) &
              == RETURN_FAIL ) goto 700
      end do
      if( write_string_tree_to_table(tree, table, err) == RETURN_FAIL ) &
           goto 800
      if( kill_string_tree(tree, err) == RETURN_FAIL ) goto 800
      ! normal exit
      table%is_null = .false.
      table%type = "string"
      answer = RETURN_SUCCESS
      return
      ! error traps
600   call err_handle(err, 1, &
           comment = "Input argument vec has length zero" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
650   call err_handle(err, 1, &
           comment = "Length of vec(1) greater than max_string_length")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
700   write(sInt,"(I12)") i
      sInt = adjustl(sInt)
      call err_handle(err, 1, &
           comment = "Failed to tabulate data" )
      call err_handle(err, 3, iobs = i )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      ijunk = nullify_table(table,err)
      return
    end function tabulate_string_variable
   !################################################################
   recursive integer(our_int) function tabulate_real_number(t, number, &
        err) result(answer)
      implicit none
      ! args
      type(dnode), pointer :: t
      real(kind=our_dble), intent(in) :: number
      type(error_type), intent(inout) :: err
      ! locals
      type(dnode), pointer :: tmp
      integer(kind=our_int) :: status
      character(len=*), parameter :: subname = "tabulate_real_number"
      ! begin
      answer = RETURN_FAIL
      tmp => t
      ! if subtree is empty, add a node
      if( .not.associated(tmp) ) then
         allocate(tmp, stat=status)
         if( status /= 0 ) goto 800
         tmp%value = number
         tmp%frequency = 1
      else if( number == tmp%value ) then
         tmp%frequency = tmp%frequency + 1
      else if( number < tmp%value ) then
         if( tabulate_real_number( tmp%left, number, err) &
              == RETURN_FAIL ) goto 999
      else
         if( tabulate_real_number( tmp%right, number, err) &
              == RETURN_FAIL ) goto 999
      end if
      t => tmp
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
999   continue
      return  
   end function tabulate_real_number
   !################################################################
   recursive integer(our_int) function tabulate_integer_number(t, &
        number, err) result(answer)
      implicit none
      ! args
      type(inode), pointer :: t
      integer(kind=our_int), intent(in) :: number
      type(error_type), intent(inout) :: err
      ! locals
      type(inode), pointer :: tmp
      integer(kind=our_int) :: status
      character(len=*), parameter :: subname = "tabulate_integer_number"
      ! begin
      answer = RETURN_FAIL
      tmp => t
      ! if subtree is empty, add a node
      if( .not.associated(tmp) ) then
         allocate(tmp, stat=status)
         if( status /= 0 ) goto 800
         tmp%value = number
         tmp%frequency = 1
      else if( number == tmp%value ) then
         tmp%frequency = tmp%frequency + 1
      else if( number < tmp%value ) then
         if( tabulate_integer_number( tmp%left, number, err) &
              == RETURN_FAIL ) goto 999
      else
         if( tabulate_integer_number( tmp%right, number, err) &
              == RETURN_FAIL ) goto 999
      end if
      t => tmp
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
999   continue
      return  
   end function tabulate_integer_number
   !################################################################
   recursive integer(our_int) function tabulate_string(t, &
        str, err) result(answer)
      implicit none
      ! args
      type(snode), pointer :: t
      character(len=*), intent(in) :: str
      type(error_type), intent(inout) :: err
      ! locals
      type(snode), pointer :: tmp
      integer(kind=our_int) :: status
      character(len=*), parameter :: subname = "tabulate_string"
      ! begin
      answer = RETURN_FAIL
      tmp => t
      ! if subtree is empty, add a node
      if( .not.associated(tmp) ) then
         allocate(tmp, stat=status)
         if( status /= 0 ) goto 800
         tmp%value = str
         tmp%frequency = 1
      else if( str == tmp%value ) then
         tmp%frequency = tmp%frequency + 1
      else if( llt ( str, tmp%value ) ) then
         if( tabulate_string( tmp%left, str, err) &
              == RETURN_FAIL ) goto 999
      else
         if( tabulate_string( tmp%right, str, err) &
              == RETURN_FAIL ) goto 999
      end if
      t => tmp
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
999   continue
      return  
   end function tabulate_string
   !################################################################
   recursive integer(our_int) function kill_real_tree(t, err) &
        result(answer)
      implicit none
      ! args
      type(dnode), pointer :: t
      type(error_type), intent(inout) :: err
      ! locals
      integer :: status
      character(len=*), parameter :: subname = "kill_real_tree"
      ! begin
      answer = RETURN_FAIL
      if(associated(t%left)) then
         if( kill_real_tree(t%left, err) == RETURN_FAIL ) goto 999
      end if
      if(associated(t%right)) then
         if( kill_real_tree(t%right, err) == RETURN_FAIL ) goto 999
      end if
      if( (.not.associated(t%left) ) .and. &
           (.not.associated(t%right)) ) then
         ! wipe out terminal node
         deallocate(t, stat=status)
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
999   continue
      return
   end function kill_real_tree
   !################################################################
   recursive integer(our_int) function kill_integer_tree(t, err) &
        result(answer)
      implicit none
      ! args
      type(inode), pointer :: t
      type(error_type), intent(inout) :: err
      ! locals
      integer :: status
      character(len=*), parameter :: subname = "kill_integer_tree"
      ! begin
      answer = RETURN_FAIL
      if(associated(t%left)) then
         if( kill_integer_tree(t%left, err) == RETURN_FAIL ) goto 999
      end if
      if(associated(t%right)) then
         if( kill_integer_tree(t%right, err) == RETURN_FAIL ) goto 999
      end if
      if( (.not.associated(t%left) ) .and. &
           (.not.associated(t%right)) ) then
         ! wipe out terminal node
         deallocate(t, stat=status)
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
999   continue
      return
   end function kill_integer_tree
   !################################################################
   recursive integer(our_int) function kill_string_tree(t, err) &
        result(answer)
      implicit none
      ! args
      type(snode), pointer :: t
      type(error_type), intent(inout) :: err
      ! locals
      integer :: status
      character(len=*), parameter :: subname = "kill_string_tree"
      ! begin
      answer = RETURN_FAIL
      if(associated(t%left)) then
         if( kill_string_tree(t%left, err) == RETURN_FAIL ) goto 999
      end if
      if(associated(t%right)) then
         if( kill_string_tree(t%right, err) == RETURN_FAIL ) goto 999
      end if
      if( (.not.associated(t%left) ) .and. &
           (.not.associated(t%right)) ) then
         ! wipe out terminal node
         deallocate(t, stat=status)
         if( status /= 0 ) goto 800
      end if
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to deallocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
999   continue
      return
    end function kill_string_tree
   !################################################################
   integer(our_int) function write_real_tree_to_table(t, table, err) &
        result(answer)
      implicit none
      ! args
      type(dnode), pointer :: t
      type(table_type), intent(inout) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "write_real_tree_to_table"
      ! go
      answer = RETURN_FAIL
      if( nullify_table( table, err) == RETURN_FAIL ) goto 800
      ! calculate tree length and allocate the arrays
      call count_dnodes(t,table)
      if( dyn_alloc( table%freq, table%length, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_alloc( table%rVal, table%length, err) &
           == RETURN_FAIL ) goto 800
      ! write values in tree to arrays
      call write_real_tree(t,table)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function write_real_tree_to_table
   !################################################################
   integer(our_int) function write_integer_tree_to_table(t, table, err) &
        result(answer)
      implicit none
      ! args
      type(inode), pointer :: t
      type(table_type), intent(inout) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "write_integer_tree_to_table"
      ! go
      answer = RETURN_FAIL
      if( nullify_table( table, err) == RETURN_FAIL ) goto 800
      ! calculate tree length and allocate the arrays
      call count_inodes(t,table)
      if( dyn_alloc( table%freq, table%length, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_alloc( table%iVal, table%length, err) &
           == RETURN_FAIL ) goto 800
      ! write values in tree to arrays
      call write_integer_tree(t,table)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function write_integer_tree_to_table
   !################################################################
   integer(our_int) function write_string_tree_to_table(t, table, err) &
        result(answer)
      implicit none
      ! args
      type(snode), pointer :: t
      type(table_type), intent(inout) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "write_string_tree_to_table"
      ! go
      answer = RETURN_FAIL
      if( nullify_table( table, err) == RETURN_FAIL ) goto 800
      ! calculate tree length and allocate the arrays
      call count_snodes(t,table)
      if( dyn_alloc( table%freq, table%length, err) &
           == RETURN_FAIL ) goto 800
      if( dyn_alloc( table%sVal, table%length, err) &
           == RETURN_FAIL ) goto 800
      ! write values in tree to arrays
      call write_string_tree(t,table)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error trap
800   continue
      call err_handle(err, 1, &
           comment = "Unable to allocate memory for object" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function write_string_tree_to_table
   !################################################################
   recursive subroutine write_real_tree(t,table)
      implicit none
      type(dnode), pointer :: t
      type(table_type), intent(inout) :: table
      if(associated(t)) then
         call write_real_tree(t%left,table)
         table%write_posn = table%write_posn + 1
         table%rVal(table%write_posn) = t%value
         table%freq(table%write_posn) = t%frequency
         call write_real_tree(t%right,table)
      end if
   end subroutine write_real_tree
   !################################################################
   recursive subroutine write_integer_tree(t,table)
      implicit none
      type(inode), pointer :: t
      type(table_type), intent(inout) :: table
      if(associated(t)) then
         call write_integer_tree(t%left,table)
         table%write_posn = table%write_posn + 1
         table%iVal(table%write_posn) = t%value
         table%freq(table%write_posn) = t%frequency
         call write_integer_tree(t%right,table)
      end if
   end subroutine write_integer_tree
   !################################################################
   recursive subroutine write_string_tree(t,table)
      implicit none
      type(snode), pointer :: t
      type(table_type), intent(inout) :: table
      if(associated(t)) then
         call write_string_tree(t%left,table)
         table%write_posn = table%write_posn + 1
         table%sVal(table%write_posn) = t%value
         table%freq(table%write_posn) = t%frequency
         call write_string_tree(t%right,table)
      end if
    end subroutine write_string_tree
   !################################################################
   recursive subroutine count_dnodes(t,table)
      implicit none
      type(dnode), pointer :: t
      type(table_type), intent(inout) :: table
      if( associated(t) ) then
         table%length = table%length + 1
         call count_dnodes(t%left,table)
         call count_dnodes(t%right,table)
      end if
   end subroutine count_dnodes
   !################################################################
   recursive subroutine count_inodes(t,table)
      implicit none
      type(inode), pointer :: t
      type(table_type), intent(inout) :: table
      if( associated(t) ) then
         table%length = table%length + 1
         call count_inodes(t%left,table)
         call count_inodes(t%right,table)
      end if
   end subroutine count_inodes
   !################################################################
   recursive subroutine count_snodes(t,table)
      implicit none
      type(snode), pointer :: t
      type(table_type), intent(inout) :: table
      if( associated(t) ) then
         table%length = table%length + 1
         call count_snodes(t%left,table)
         call count_snodes(t%right,table)
      end if
    end subroutine count_snodes
   !################################################################
   integer(our_int) function get_table_length(length, table, err) &
        result(answer)
      implicit none
      ! args
      integer(our_int), intent(out) :: length
      type(table_type), intent(in) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_table_length"
      ! go
      answer = RETURN_FAIL
      if( table%is_null ) goto 700
      length = table%length
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Table object is null")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function get_table_length
   !################################################################
   integer(our_int) function get_table_type(ttype, table, err) &
        result(answer)
      implicit none
      ! args
      character(len=*), intent(out) :: ttype
      type(table_type), intent(in) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_table_type"
      ! go
      answer = RETURN_FAIL
      ttype = ""
      if( table%is_null ) goto 700
      ttype = table%type
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Table object is null")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function get_table_type
   !################################################################
   integer(our_int) function get_real_table_values(val, table, err) &
        result(answer)
      implicit none
      ! args
      real(kind=our_dble), pointer :: val(:)
      type(table_type), intent(in) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_real_table_values"
      ! go
      answer = RETURN_FAIL
      if( table%is_null ) goto 700
      if( table%type /= "double precision" ) goto 750
      if( dyn_alloc(val, table%length, err) == RETURN_FAIL ) goto 800
      val(:) = table%rVal(:)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Table object is null")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
750   call err_handle(err, 1, &
           comment = "Output argument does not match table type" )
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function get_real_table_values
   !################################################################
   integer(our_int) function get_integer_table_values(val, table, err) &
        result(answer)
      implicit none
      ! args
      integer(kind=our_int), pointer :: val(:)
      type(table_type), intent(in) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_integer_table_values"
      ! go
      answer = RETURN_FAIL
      if( table%is_null ) goto 700
      if( table%type /= "integer" ) goto 750
      if( dyn_alloc(val, table%length, err) == RETURN_FAIL ) goto 800
      val(:) = table%iVal(:)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Table object is null")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
750   call err_handle(err, 1, &
           comment = "Output argument does not match table type")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function get_integer_table_values
   !################################################################
   integer(our_int) function get_string_table_values(val, table, err) &
        result(answer)
      implicit none
      ! args
      character(len=*), pointer :: val(:)
      type(table_type), intent(in) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_string_table_values"
      ! go
      answer = RETURN_FAIL
      if( table%is_null ) goto 700
      if( table%type /= "string" ) goto 750
      if( dyn_alloc(val, table%length, err) == RETURN_FAIL ) goto 800
      val(:) = table%sVal(:)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Table object is null")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
750   call err_handle(err, 1, &
           comment = "Output argument does not match table type")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
    end function get_string_table_values
   !################################################################
   integer(our_int) function get_table_frequencies(freq, table, err) &
        result(answer)
      implicit none
      ! args
      integer(kind=our_int), pointer :: freq(:)
      type(table_type), intent(in) :: table
      type(error_type), intent(inout) :: err
      ! locals
      character(len=*), parameter :: subname = "get_table_frequencies"
      ! go
      answer = RETURN_FAIL
      if( table%is_null ) goto 700
      if( dyn_alloc(freq, table%length, err) == RETURN_FAIL ) goto 800
      freq(:) = table%freq(:)
      ! normal exit
      answer = RETURN_SUCCESS
      return
      ! error traps
700   call err_handle(err, 1, &
           comment = "Table object is null")
      call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
800   call err_handle(err, 2, whichsub = subname, whichmod = modname )
      return
   end function get_table_frequencies
   !################################################################
end module tabulate
!#######################################################################
