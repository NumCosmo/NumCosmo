AC_DEFUN([AX_GFORTRAN_VERSION], [
  AC_LANG_PUSH([Fortran])
  GFORTRAN_VERSION=""

  AS_IF([test "x$ac_cv_fc_compiler_gnu" = "xyes"],[
    AS_IF([$FC -dumpversion > /dev/null 2>&1],[
      AC_CACHE_CHECK([gfortran version],[ax_cv_gfortran_version],[
        ax_cv_gfortran_version="`$FC -dumpversion`"
        AS_IF([test "x$ax_cv_gfortran_version" = "x"],[
          ax_cv_gfortran_version=""
        ])
      ])
      GFORTRAN_VERSION=$ax_cv_gfortran_version
    ])
  ])
  AC_SUBST([GFORTRAN_VERSION])
  AC_LANG_POP([Fortran])
])
