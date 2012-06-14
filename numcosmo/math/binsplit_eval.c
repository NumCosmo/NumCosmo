
#ifndef NC_BINSPLIT_EVAL_NAME
#error "To include binsplit_eval.c you must define NC_BINSPLIT_EVAL_NAME"
#endif

#if !defined(_BINSPLIT_FUNC_P) || !defined(_BINSPLIT_FUNC_Q)
#error "To include binsplit_eval.c you must define at least the macros _BINSPLIT_FUNC_(P|Q)"
#endif

void
NC_BINSPLIT_EVAL_NAME (NcmBinSplit *bs, gulong n1, gulong n2)
{
  gulong nd = n2 - n1;
  bs->n1 = n1;
  bs->n2 = n2;
  //printf ("# nd [%lu %lu) %lu %d\n", n1, n2, nd, nd < 5);
  if (nd < 5)
  {
    _BINSPLIT_FUNC_P (bs->P, NCM_BINSPLIT_ONE, n1, bs->userdata);
    _BINSPLIT_FUNC_Q (bs->Q, NCM_BINSPLIT_ONE, n1, bs->userdata);
    _BINSPLIT_FUNC_B (bs->B, NCM_BINSPLIT_ONE, n1, bs->userdata);

    switch (nd)
    {
      case 1:
        mpz_set (bs->T, bs->P);
        _BINSPLIT_FUNC_A (bs->T, bs->T, n1, bs->userdata);
        break;
      case 2:
        _BINSPLIT_FUNC_P (bs->P, bs->P, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->Q, bs->Q, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->B, bs->B, n1+1L, bs->userdata);

        _BINSPLIT_FUNC_P (bs->temp1, NCM_BINSPLIT_ONE, n1, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp1, bs->temp1, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_A (bs->temp1, bs->temp1, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp1, bs->temp1, n1+1L, bs->userdata);

        mpz_set (bs->temp2, bs->P);
        _BINSPLIT_FUNC_A (bs->temp2, bs->temp2, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp2, bs->temp2, n1, bs->userdata);

        mpz_add (bs->T, bs->temp1, bs->temp2);
        break;
      case 3:
        _BINSPLIT_FUNC_P (bs->P, bs->P, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_P (bs->P, bs->P, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->Q, bs->Q, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->Q, bs->Q, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->B, bs->B, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->B, bs->B, n1+2L, bs->userdata);

        _BINSPLIT_FUNC_P (bs->temp1, NCM_BINSPLIT_ONE, n1, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp1, bs->temp1, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp1, bs->temp1, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_A (bs->temp1, bs->temp1, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp1, bs->temp1, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp1, bs->temp1, n1+2L, bs->userdata);

        _BINSPLIT_FUNC_P (bs->temp2, NCM_BINSPLIT_ONE, n1, bs->userdata);
        _BINSPLIT_FUNC_P (bs->temp2, bs->temp2, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp2, bs->temp2, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_A (bs->temp2, bs->temp2, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp2, bs->temp2, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp2, bs->temp2, n1+2L, bs->userdata);

        mpz_set (bs->temp3, bs->P);
        _BINSPLIT_FUNC_A (bs->temp3, bs->temp3, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp3, bs->temp3, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp3, bs->temp3, n1+1L, bs->userdata);

        mpz_add (bs->T, bs->temp1, bs->temp2);
        mpz_add (bs->T, bs->T, bs->temp3);
        break;
      case 4:
        _BINSPLIT_FUNC_P (bs->P, bs->P, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_P (bs->P, bs->P, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_P (bs->P, bs->P, n1+3L, bs->userdata);        
        _BINSPLIT_FUNC_Q (bs->Q, bs->Q, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->Q, bs->Q, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->Q, bs->Q, n1+3L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->B, bs->B, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->B, bs->B, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->B, bs->B, n1+3L, bs->userdata);

        _BINSPLIT_FUNC_P (bs->temp1, NCM_BINSPLIT_ONE, n1, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp1, bs->temp1, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp1, bs->temp1, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp1, bs->temp1, n1+3L, bs->userdata);
        _BINSPLIT_FUNC_A (bs->temp1, bs->temp1, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp1, bs->temp1, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp1, bs->temp1, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp1, bs->temp1, n1+3L, bs->userdata);

        _BINSPLIT_FUNC_P (bs->temp2, NCM_BINSPLIT_ONE, n1, bs->userdata);
        _BINSPLIT_FUNC_P (bs->temp2, bs->temp2, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp2, bs->temp2, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp2, bs->temp2, n1+3L, bs->userdata);
        _BINSPLIT_FUNC_A (bs->temp2, bs->temp2, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp2, bs->temp2, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp2, bs->temp2, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp2, bs->temp2, n1+3L, bs->userdata);

        _BINSPLIT_FUNC_P (bs->temp3, NCM_BINSPLIT_ONE, n1, bs->userdata);
        _BINSPLIT_FUNC_P (bs->temp3, bs->temp3, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_P (bs->temp3, bs->temp3, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_Q (bs->temp3, bs->temp3, n1+3L, bs->userdata);
        _BINSPLIT_FUNC_A (bs->temp3, bs->temp3, n1+2L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp3, bs->temp3, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp3, bs->temp3, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp3, bs->temp3, n1+3L, bs->userdata);
        
        mpz_set (bs->temp4, bs->P);
        _BINSPLIT_FUNC_A (bs->temp4, bs->temp4, n1+3L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp4, bs->temp4, n1, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp4, bs->temp4, n1+1L, bs->userdata);
        _BINSPLIT_FUNC_B (bs->temp4, bs->temp4, n1+2L, bs->userdata);

        mpz_add (bs->T, bs->temp1, bs->temp2);
        mpz_add (bs->T, bs->T, bs->temp3);
        mpz_add (bs->T, bs->T, bs->temp4);
        break;
      default:
        g_error ("binsplit_eval[%s] should not be called with (%ld, %ld) nd = %ld\n", __func__, n1, n2, nd);
        break;
    }
  }
  else
  {
    gulong nl = (n2 + n1) / 2;

    if (bs->bs[0] == NULL)
      bs->bs[0] = ncm_binsplit_alloc (bs->userdata);
    if (bs->bs[1] == NULL)
      bs->bs[1] = ncm_binsplit_alloc (bs->userdata);

    NC_BINSPLIT_EVAL_NAME (bs->bs[0], n1, nl);
    NC_BINSPLIT_EVAL_NAME (bs->bs[1], nl, n2);

    mpz_mul (bs->P, bs->bs[0]->P, bs->bs[1]->P);
    mpz_mul (bs->Q, bs->bs[0]->Q, bs->bs[1]->Q);
#ifdef _HAS_FUNC_B
    mpz_mul (bs->B, bs->bs[0]->B, bs->bs[1]->B);
#endif
    mpz_mul (bs->temp1, bs->bs[0]->P, bs->bs[1]->T);
    mpz_mul (bs->temp2, bs->bs[0]->T, bs->bs[1]->Q);

#ifdef _HAS_FUNC_B
    mpz_mul (bs->temp1, bs->temp1, bs->bs[0]->B);    
    mpz_mul (bs->temp2, bs->temp2, bs->bs[1]->B);
#endif
    
    mpz_add (bs->T, bs->temp1, bs->temp2);
  }
//  mpfr_printf ("# NUC %.15g %Zd/(%Zd %Zd)\n", ncm_binsplit_get_d (bs, GMP_RNDD), bs->T, bs->Q, bs->B);
}

#undef NC_BINSPLIT_EVAL_NAME
#undef _BINSPLIT_FUNC_P
#undef _BINSPLIT_FUNC_Q
#undef _BINSPLIT_FUNC_A
#undef _BINSPLIT_FUNC_B
#undef _HAS_FUNC_B
