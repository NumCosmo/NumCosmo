! ============================================================================
! "Read_Archive_Map" reads a sky map from a binary table in a FITS file.  The
! map consists of two columns in the table; the first contains temperature
! and the second contains N_Obs.
!
! Arguments:
!   FileName  Character  The FITS filename.
!   Temp      Real*4     The array of temperatures.
!   NObs      Real*4     The array of observation numbers.
!   NPix      Integer    The number of elements in the output arrays.
!   Status    Integer    A status code: 0=success.
!
! Written by Michael R. Greason, SSAI, 07 February 2006.
! ============================================================================
Subroutine Read_Archive_Map (FileName, Temp, NObs, NPix, Status)
!
Character(*),                Intent(In)  :: FileName
Integer (Kind=4),            Intent(Out) :: NPix
Integer (Kind=4),            Intent(Out) :: Status
Real (Kind=4), Dimension(:), Intent(Out) :: Temp
Real (Kind=4), Dimension(:), Intent(Out) :: NObs
!
Character (80)   :: str
Integer (Kind=4) :: i, n, Unit, blk, typ
Logical          :: anyf
! ----------------------------------------------------------------------------
!           Open the FITS file.
!
NPix   = 0
Status = 0
i      = 0
Call FTGIOU (Unit, Status)
Call FTOPEN (Unit, FileName, i, blk, Status)
!
!           Select the HDU.
!
i   = 1
typ = 0
Do While ((typ .NE. 2) .AND. (Status .EQ. 0))
  i = i + 1
  Call FTMAHD (Unit, i, typ, Status)
End Do
If ((typ .NE. 2) .OR. (Status .NE. 0)) Then
  Print 1, trim(FileName)
  Return
End If
!
!           Extract the number of pixels in the map.
!
Call FTGNRW (Unit, i, Status)
If (Status .NE. 0) Return
NPix = i
!
!           Extract the map.
!
Call FTGCVE (Unit, 1, 1, 1, NPix, 0.0E0, Temp, anyf, Status)
Call FTGCVE (Unit, 2, 1, 1, NPix, 0.0E0, NObs, anyf, Status)
!
!           Close the FITS file.
!
Call FTCLOS (Unit, Status)
If (Unit .GE. 50) Call FTFIOU (Unit, Status)
!
Return
! ----------------------------------------------------------------------------
1 Format ('Read_Archive_Map: There is no binary FITS table extension in ',A)
!
End Subroutine Read_Archive_Map
