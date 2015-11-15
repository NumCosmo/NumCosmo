! ============================================================================
! 'Read_FITS_Real_2D' reads a two-dimensional floating point array from the
! primary header/data unit in a FITS file.  This routine allocates the array
! from the heap if it hasn't already been associated.
!
! If the array has already been created, this routine assumes that it has
! the correct size and dimensions!  If allocated here then it will be indexed
! like a C array:  0:n-1.
!
! Arguments:
!	FileName - The name of the FITS file.
!	Arr      - The array to store the data in.
!	Status   - A status code: 0=success.
!	Dim1     - The first dimension size.  Optional.
!	Dim2     - The second dimension size.  Optional.
!	IndFmt   - A bit-coded optional flag indicating how the array should
!	           be indexed in each axis; 0=0-based, 1=1-based.  Defaults
!	           to 0---all 0-based.  This is only relevant if this routine
!	           allocates the array.
!
! Written by Michael R. Greason, SSAI, 07 February 2006.
! ============================================================================
Subroutine Read_FITS_Real_2D (FileName, Arr, Status, Dim1, Dim2, IndFmt)
!
Character(*),               Intent(In)  :: FileName
Real (Kind=4),  Pointer, Dimension(:,:) :: Arr
Integer (Kind=4),           Intent(Out) :: Status
Integer (Kind=4), Optional, Intent(Out) :: Dim1
Integer (Kind=4), Optional, Intent(Out) :: Dim2
Integer (Kind=4), Optional, Intent(In)  :: IndFmt
!
Character (80)   :: Comm
Integer (Kind=4) :: Unit, Blk, I, D1, D2, NEle, Ifmt
Logical          :: AnyF
! ----------------------------------------------------------------------------
Ifmt = 0
If (Present(IndFmt)) Ifmt = IndFmt
!
!			Open the FITS file.
!
Status = 0
I      = 0
Call FTGIOU (Unit, Status)
Call FTOPEN (Unit, FileName, I, Blk, Status)
!
!			Allocate space for the array in the primary HDU.
!
Call FTGIDM (Unit, D1, Status)
If ((Status .NE. 0) .OR. (D1 .LT. 2)) Then
  If (D1 .LT. 2) Print 1, Status, 'Primary HDU has insufficient axes.'
  Return
End If
!
Call FTGKYJ (Unit, 'NAXIS1', D1, Comm, Status)
Call FTGKYJ (Unit, 'NAXIS2', D2, Comm, Status)
If (Status .NE. 0) Return
If (Present(Dim1)) Dim1 = D1
If (Present(Dim2)) Dim2 = D2
!
If (.NOT. Associated(Arr)) Then
  Select Case (Ifmt)
    Case (3)
      Allocate (Arr(1:D1,     1:D2    ), Stat=Status)
    Case (2)
      Allocate (Arr(0:(D1-1), 1:D2    ), Stat=Status)
    Case (1)
      Allocate (Arr(1:D1,     0:(D2-1)), Stat=Status)
    Case Default
      Allocate (Arr(0:(D1-1), 0:(D2-1)), Stat=Status)
  End Select
  If (Status .NE. 0) Return
End If
!
!			Fill the array.
!
NELE = D1 * D2
Call FTGPVE (Unit, 0, 1, NEle, 0, Arr, AnyF, Status)
!
!			Close the FITS file.
!
Call FTCLOS (Unit, Status)
If (Unit .GE. 50) Call FTFIOU (Unit, Status)
!
Return
! ----------------------------------------------------------------------------
1 Format ('Read_FITS_Real_2D:  Status = ',I10,3X,A)
!
End Subroutine Read_FITS_Real_2D
! ============================================================================
! 'Read_FITS_Real_2D' reads a two-dimensional floating point array from the
! primary header/data unit in a FITS file.  This routine allocates the array
! from the heap if it hasn't already been associated.
!
! If the array has already been created, this routine assumes that it has
! the correct size and dimensions!  If allocated here then it will be indexed
! like a C array:  0:n-1.
!
! Arguments:
!	FileName - The name of the FITS file.
!	Arr      - The array to store the data in.
!	Status   - A status code: 0=success.
!	Dim1     - The first dimension size.  Optional.
!	Dim2     - The second dimension size.  Optional.
!	IndFmt   - A bit-coded optional flag indicating how the array should
!	           be indexed in each axis; 0=0-based, 1=1-based.  Defaults
!	           to 0---all 0-based.  This is only relevant if this routine
!	           allocates the array.
!
! Written by Michael R. Greason, SSAI, 07 February 2006.
! ============================================================================
Subroutine Read_FITS_Double_2D (FileName, Arr, Status, Dim1, Dim2, IndFmt)
!
Character(*),               Intent(In)  :: FileName
Real (Kind=8),  Pointer, Dimension(:,:) :: Arr
Integer (Kind=4),           Intent(Out) :: Status
Integer (Kind=4), Optional, Intent(Out) :: Dim1
Integer (Kind=4), Optional, Intent(Out) :: Dim2
Integer (Kind=4), Optional, Intent(In)  :: IndFmt
!
Character (80)   :: Comm
Integer (Kind=4) :: Unit, Blk, I, D1, D2, NEle, Ifmt
Logical          :: AnyF
! ----------------------------------------------------------------------------
Ifmt = 0
If (Present(IndFmt)) Ifmt = IndFmt
!
!			Open the FITS file.
!
Status = 0
I      = 0
Call FTGIOU (Unit, Status)
Call FTOPEN (Unit, FileName, I, Blk, Status)
!
!			Allocate space for the array in the primary HDU.
!
Call FTGIDM (Unit, D1, Status)
If ((Status .NE. 0) .OR. (D1 .LT. 2)) Then
  If (D1 .LT. 2) Print 1, Status, 'Primary HDU has insufficient axes.'
  Return
End If
!
Call FTGKYJ (Unit, 'NAXIS1', D1, Comm, Status)
Call FTGKYJ (Unit, 'NAXIS2', D2, Comm, Status)
If (Status .NE. 0) Return
If (Present(Dim1)) Dim1 = D1
If (Present(Dim2)) Dim2 = D2
!
If (.NOT. Associated(Arr)) Then
  Select Case (Ifmt)
    Case (3)
      Allocate (Arr(1:D1,     1:D2    ), Stat=Status)
    Case (2)
      Allocate (Arr(0:(D1-1), 1:D2    ), Stat=Status)
    Case (1)
      Allocate (Arr(1:D1,     0:(D2-1)), Stat=Status)
    Case Default
      Allocate (Arr(0:(D1-1), 0:(D2-1)), Stat=Status)
  End Select
  If (Status .NE. 0) Return
End If
!
!			Fill the array.
!
NELE = D1 * D2
Call FTGPVD (Unit, 0, 1, NEle, 0, Arr, AnyF, Status)
!
!			Close the FITS file.
!
Call FTCLOS (Unit, Status)
If (Unit .GE. 50) Call FTFIOU (Unit, Status)
!
Return
! ----------------------------------------------------------------------------
1 Format ('Read_FITS_Double_2D:  Status = ',I10,3X,A)
!
End Subroutine Read_FITS_Double_2D
! ============================================================================
! 'Read_FITS_Real_3D' reads a three-dimensional floating point array from the
! primary header/data unit in a FITS file.  This routine allocates the array
! from the heap if it hasn't already been associated.
!
! If the array has already been created, this routine assumes that it has
! the correct size and dimensions!  If allocated here then it will be indexed
! like a C array:  0:n-1.
!
! Arguments:
!	FileName - The name of the FITS file.
!	Arr      - The array to store the data in.
!	Status   - A status code: 0=success.
!	Dim1     - The first dimension size.  Optional.
!	Dim2     - The second dimension size.  Optional.
!	Dim3     - The third dimension size.  Optional.
!	IndFmt   - A bit-coded optional flag indicating how the array should
!	           be indexed in each axis; 0=0-based, 1=1-based.  Defaults
!	           to 0---all 0-based.  This is only relevant if this routine
!	           allocates the array.
!
! Written by Michael R. Greason, SSAI, 07 February 2006.
! ============================================================================
Subroutine Read_FITS_Real_3D (FileName, Arr, Status, Dim1, Dim2, Dim3, IndFmt)
!
Character(*),                Intent(In)  :: FileName
Real (Kind=4), Pointer, Dimension(:,:,:) :: Arr
Integer (Kind=4),            Intent(Out) :: Status
Integer (Kind=4),  Optional, Intent(Out) :: Dim1
Integer (Kind=4),  Optional, Intent(Out) :: Dim2
Integer (Kind=4),  Optional, Intent(Out) :: Dim3
Integer (Kind=4),  Optional, Intent(In)  :: IndFmt
!
Character (80)   :: Comm
Integer (Kind=4) :: Unit, Blk, I, D1, D2, D3, NEle, Ifmt
Logical          :: AnyF
! ----------------------------------------------------------------------------
Ifmt = 0
If (Present(IndFmt)) Ifmt = IndFmt
!
!			Open the FITS file.
!
Status = 0
I      = 0
Call FTGIOU (Unit, Status)
Call FTOPEN (Unit, FileName, I, Blk, Status)
!
!			Allocate space for the array in the primary HDU.
!
Call FTGIDM (Unit, D1, Status)
If ((Status .NE. 0) .OR. (D1 .LT. 3)) Then
  If (D1 .LT. 3) Print 1, Status, 'Primary HDU has insufficient axes.'
  Return
End If
!
Call FTGKYJ (Unit, 'NAXIS1', D1, Comm, Status)
Call FTGKYJ (Unit, 'NAXIS2', D2, Comm, Status)
Call FTGKYJ (Unit, 'NAXIS3', D3, Comm, Status)
If (Status .NE. 0) Return
If (Present(Dim1)) Dim1 = D1
If (Present(Dim2)) Dim2 = D2
If (Present(Dim3)) Dim3 = D3
!
If (.NOT. Associated(Arr)) Then
  Select Case (Ifmt)
    Case (7)
      Allocate (Arr(1:D1,     1:D2,     1:D3    ), Stat=Status)
    Case (6)
      Allocate (Arr(0:(D1-1), 1:D2,     1:D3    ), Stat=Status)
    Case (5)
      Allocate (Arr(1:D1,     0:(D2-1), 1:D3    ), Stat=Status)
    Case (4)
      Allocate (Arr(0:(D1-1), 0:(D2-1), 1:D3    ), Stat=Status)
    Case (3)
      Allocate (Arr(1:D1,     1:D2,     0:(D3-1)), Stat=Status)
    Case (2)
      Allocate (Arr(0:(D1-1), 1:D2,     0:(D3-1)), Stat=Status)
    Case (1)
      Allocate (Arr(1:D1,     0:(D2-1), 0:(D3-1)), Stat=Status)
    Case Default
      Allocate (Arr(0:(D1-1), 0:(D2-1), 0:(D3-1)), Stat=Status)
  End Select
  If (Status .NE. 0) Return
End If
!
!			Fill the array.
!
NELE = D1 * D2 * D3
Call FTGPVE (Unit, 0, 1, NEle, 0, Arr, AnyF, Status)
!
!			Close the FITS file.
!
Call FTCLOS (Unit, Status)
If (Unit .GE. 50) Call FTFIOU (Unit, Status)
!
Return
! ----------------------------------------------------------------------------
1 Format ('Read_FITS_Real_3D:  Status = ',I10,3X,A)
!
End Subroutine Read_FITS_Real_3D
! ============================================================================
! 'Read_FITS_Complex_2D' reads a two-dimensional complex array from the primary
! header/data unit in a FITS file.  This routine allocates the array from the
! heap if it hasn't already been associated.
!
! If the array has already been created, this routine assumes that it has
! the correct size and dimensions!  If allocated here then it will be indexed
! like a C array:  0:n-1.
!
! FITS does not directly support complex data; it treats each element as a
! two-element array.  Therefore, if the complex array has dimensions MxN it is
! stored as a 2xMxN array in the FITS file.  This routine will use
! 'Read_FITS_Real_3D' to read the file into a 3D array; this array will be
! used to fill the output array.
!
! Arguments:
!	FileName - The name of the FITS file.
!	Arr      - The array to store the data in.
!	Status   - A status code: 0=success.
!	Dim1     - The first dimension size.  Optional.
!	Dim2     - The second dimension size.  Optional.
!	IndFmt   - A bit-coded optional flag indicating how the array should
!	           be indexed in each axis; 0=0-based, 1=1-based.  Defaults
!	           to 0---all 0-based.
!
! Written by Michael R. Greason, SSAI, 07 February 2006.
! ============================================================================
Subroutine Read_FITS_Complex_2D (FileName, Arr, Status, Dim1, Dim2, IndFmt)
!
Character(*),                 Intent(In)  :: FileName
Complex (Kind=4), Pointer, Dimension(:,:) :: Arr
Integer (Kind=4),             Intent(Out) :: Status
Integer (Kind=4),   Optional, Intent(Out) :: Dim1
Integer (Kind=4),   Optional, Intent(Out) :: Dim2
Integer (Kind=4),   Optional, Intent(In)  :: IndFmt
!
Integer (Kind=4) :: DC, D1, D2, I, I0, I1, J, J0, J1, Ifmt
Real (Kind=4), Pointer, Dimension(:,:,:) :: Tmp
!
Interface
Subroutine Read_FITS_Real_3D (FileName, Arr, Status, Dim1, Dim2, Dim3, IndFmt)
Character(*),                Intent(In)  :: FileName
Real (Kind=4), Pointer, Dimension(:,:,:) :: Arr
Integer (Kind=4),            Intent(Out) :: Status
Integer (Kind=4),  Optional, Intent(Out) :: Dim1
Integer (Kind=4),  Optional, Intent(Out) :: Dim2
Integer (Kind=4),  Optional, Intent(Out) :: Dim3
Integer (Kind=4),  Optional, Intent(In)  :: IndFmt
End Subroutine Read_FITS_Real_3D
End Interface
! ----------------------------------------------------------------------------
Ifmt = 0
If (Present(IndFmt)) Ifmt = IndFmt
!
!			Read the file into a temporary array.
!
I = Ifmt * 2 + 1
Call Read_FITS_Real_3D (FileName, Tmp, Status, DC, D1, D2, I)
If (Status .NE. 0) Return
If ((DC .NE. 2) .OR. (D1 .LE. 0) .OR. (D2 .LE. 0)) Then
  Status = 1
  Return
End If
If (Present(Dim1)) Dim1 = D1
If (Present(Dim2)) Dim2 = D2
!
!			Allocate and fill the output array.
!
Select Case (Ifmt)
  Case (3)
    I0 = 1
    I1 = D1
    J0 = 1
    J1 = D2
  Case (2)
    I0 = 0
    I1 = D1 - 1
    J0 = 1
    J1 = D2
  Case (1)
    I0 = 1
    I1 = D1
    J0 = 0
    J1 = D2 -1
  Case Default
    I0 = 0
    I1 = D1 - 1
    J0 = 0
    J1 = D2 -1
End Select
If (.NOT. Associated(Arr)) Then
  Allocate (Arr(I0:I1,J0:J1), Stat=Status)
  If (Status .NE. 0) Return
End If
!
Do J = J0, J1
  Do I = I0, I1
    Arr(I,J) = cmplx(Tmp(1,I,J), Tmp(2,I,J), Kind=4)
  End Do
End Do
!
Deallocate(Tmp)
!
Return
! ----------------------------------------------------------------------------
End Subroutine Read_FITS_Complex_2D
! ============================================================================
! 'Read_FITS_Complex_2D_LM' reads a two-dimensional complex array from the
! primary header/data unit in a FITS file.  This routine reads and fills the
! array one element at a time in an attempt to conserve memory at the cost of
! performance.
!
! If the array has already been created, this routine assumes that it has
! the correct size and dimensions!  If allocated here then it will be indexed
! like a C array:  0:n-1.
!
! FITS does not directly support complex data; it treats each element as a
! two-element array.  Therefore, if the complex array has dimensions MxN it is
! stored as a 2xMxN array in the FITS file.
!
! Arguments:
!	FileName - The name of the FITS file.
!	Arr      - The array to store the data in.
!	Status   - A status code: 0=success.
!	Dim1     - The first dimension size.  Optional.
!	Dim2     - The second dimension size.  Optional.
!	IndFmt   - A bit-coded optional flag indicating how the array should
!	           be indexed in each axis; 0=0-based, 1=1-based.  Defaults
!	           to 0---all 0-based.
!
! Written by Michael R. Greason, SSAI, 07 February 2006.
! ============================================================================
Subroutine Read_FITS_Complex_2D_LM (FileName, Arr, Status, Dim1, Dim2, IndFmt)
!
Character(*),                 Intent(In)  :: FileName
Complex (Kind=4), Pointer, Dimension(:,:) :: Arr
Integer (Kind=4),             Intent(Out) :: Status
Integer (Kind=4),   Optional, Intent(Out) :: Dim1
Integer (Kind=4),   Optional, Intent(Out) :: Dim2
Integer (Kind=4),   Optional, Intent(In)  :: IndFmt
!
Character (80)   :: Comm
Integer (Kind=4) :: D1, D2, I, I0, I1, J, J0, J1, K!, Ifmt
Integer (Kind=4) :: Unit, Blk, Ifmt
Logical          :: AnyF
Real (Kind=4), Dimension(2) :: Tmp
! ----------------------------------------------------------------------------
Ifmt = 0
If (Present(IndFmt)) Ifmt = IndFmt
!
!			Open the FITS file.
!
Status = 0
I      = 0
Call FTGIOU (Unit, Status)
Call FTOPEN (Unit, FileName, I, Blk, Status)
!
!			Allocate space for the array in the primary HDU.
!
Call FTGIDM (Unit, D1, Status)
If ((Status .NE. 0) .OR. (D1 .LT. 3)) Then
  If (D1 .LT. 3) Print 1, Status, 'Primary HDU has insufficient axes.'
  Return
End If
!
Call FTGKYJ (Unit, 'NAXIS2', D1, Comm, Status)
Call FTGKYJ (Unit, 'NAXIS3', D2, Comm, Status)
If (Status .NE. 0) Return
If (Present(Dim1)) Dim1 = D1
If (Present(Dim2)) Dim2 = D2
!
Select Case (Ifmt)
  Case (3)
    I0 = 1
    I1 = D1
    J0 = 1
    J1 = D2
  Case (2)
    I0 = 0
    I1 = D1 - 1
    J0 = 1
    J1 = D2
  Case (1)
    I0 = 1
    I1 = D1
    J0 = 0
    J1 = D2 - 1
  Case Default
    I0 = 0
    I1 = D1 - 1
    J0 = 0
    J1 = D2 - 1
End Select
If (.NOT. Associated(Arr)) Then
  Allocate (Arr(I0:I1,J0:J1), Stat=Status)
  If (Status .NE. 0) Return
End If
!
!			Fill the array.
!
K = -1
CmplxLoop:							&
Do J = J0, J1
  Do I = I0, I1
    K = K + 2
    Call FTGPVE (Unit, 0, K, 2, 0, Tmp, AnyF, Status)
    If (Status .NE. 0) Exit CmplxLoop
    Arr(I,J) = cmplx(Tmp(1), Tmp(2), Kind=4)
  End Do
End Do CmplxLoop
!
!			Close the FITS file.
!
Call FTCLOS (Unit, Status)
If (Unit .GE. 50) Call FTFIOU (Unit, Status)
!
Return
! ----------------------------------------------------------------------------
1 Format ('Read_FITS_Complex_2D_LM:  Status = ',I10,3X,A)
!
Return
! ----------------------------------------------------------------------------
End Subroutine Read_FITS_Complex_2D_LM
