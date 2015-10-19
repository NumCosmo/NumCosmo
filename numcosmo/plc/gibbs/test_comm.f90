program test_comm

	use comm_br_mod
	integer::handle
	call comm_br_initialize_object('/users/benabed/Desktop/sigma_dx9_delta_WMAP_n128_v1.fits','/users/benabed/Desktop/test_cls_v4.dat', 2, 32, 1, 20, 20, 100, 1, handle)


end
