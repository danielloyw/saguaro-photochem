  INTEGER FUNCTION FIND_NAME(xname,name)
    USE PRECISION
    IMPLICIT none
    CHARACTER(len=12), INTENT(IN) :: xname
    CHARACTER(len=12), INTENT(IN), DIMENSION(0:) :: name
    LOGICAL found
    INTEGER :: m, nsp
    nsp = SIZE(name,1)-1
    found = .false.
    m = 0
    DO WHILE ((.not. found).and.(m <= nsp))
       IF(TRIM(ADJUSTL(xname)) == TRIM(ADJUSTL(name(m)))) THEN
          FIND_NAME = m
          found = .true.
          EXIT
       END IF
       m = m + 1
    END DO
    IF(.not. found) THEN
       FIND_NAME = 0
    END IF
    RETURN
  END FUNCTION FIND_NAME
