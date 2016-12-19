DOUBLE PRECISION FUNCTION DGAMMA(X)
      INTEGER I,N
      LOGICAL PARITY
      DIMENSION C(7),P(8),Q(8)
      DOUBLE PRECISION                                           &
          C,CONV,EPS,FACT,HALF,ONE,P,PI,Q,RES,SQRTPI,SUM,TWELVE, &
          TWO,X,XBIG,XDEN,XINF,XMININ,XNUM,Y,Y1,YSQ,Z,ZERO  ! gamma.irp.f:  95
      DATA ONE,HALF,TWELVE,TWO,ZERO/1.0D0,0.5D0,12.0D0,2.0D0,0.0D0/,&
           SQRTPI/0.9189385332046727417803297D0/, &
           PI/3.1415926535897932384626434D0/  ! gamma.irp.f: 105
      DATA XBIG,XMININ,EPS/171.624D0,2.23D-308,2.22D-16/,&
           XINF/1.79D308/  ! gamma.irp.f: 113
      DATA P/-1.71618513886549492533811D+0,2.47656508055759199108314D+1,&
             -3.79804256470945635097577D+2,6.29331155312818442661052D+2, &
             8.66966202790413211295064D+2,-3.14512729688483675254357D+4, &
             -3.61444134186911729807069D+4,6.64561438202405440627855D+4/  ! gamma.irp.f: 127
      DATA Q/-3.08402300119738975254353D+1,3.15350626979604161529144D+2,&
            -1.01515636749021914166146D+3,-3.10777167157231109440444D+3, &
              2.25381184209801510330112D+4,4.75584627752788110767815D+3, &
            -1.34659959864969306392456D+5,-1.15132259675553483497211D+5/  ! gamma.irp.f: 131
      DATA C/-1.910444077728D-03,8.4171387781295D-04,&
           -5.952379913043012D-04,7.93650793500350248D-04, &
           -2.777777777777681622553D-03,8.333333333333333331554247D-02, &
            5.7083835261D-03/  ! gamma.irp.f: 142
           I  = DBLE(I)                                               ! gamma.irp.f: 150
      PARITY = .FALSE.                                                ! gamma.irp.f: 151
      FACT = ONE                                                      ! gamma.irp.f: 152
      N = 0                                                           ! gamma.irp.f: 153
      Y = X                                                           ! gamma.irp.f: 154
      IF (Y .LE. ZERO) THEN                                           ! gamma.irp.f: 155
            Y = -X                                                    ! gamma.irp.f: 159
            Y1 = AINT(Y)                                              ! gamma.irp.f: 160
            RES = Y - Y1                                              ! gamma.irp.f: 161
            IF (RES .NE. ZERO) THEN                                   ! gamma.irp.f: 162
                  IF (Y1 .NE. AINT(Y1*HALF)*TWO) then                 ! gamma.irp.f: 163
 PARITY = .TRUE.                                                      ! gamma.irp.f: 163
  endif                                                               ! gamma.irp.f: 163
                  FACT = -PI / SIN(PI*RES)                            ! gamma.irp.f: 164
                  Y = Y + ONE                                         ! gamma.irp.f: 165
               ELSE                                                   ! gamma.irp.f: 166
                  RES = XINF                                          ! gamma.irp.f: 167
                  GO TO 900                                           ! gamma.irp.f: 168
            END IF                                                    ! gamma.irp.f: 169
      END IF                                                          ! gamma.irp.f: 170
      IF (Y .LT. EPS) THEN                                            ! gamma.irp.f: 174
            IF (Y .GE. XMININ) THEN                                   ! gamma.irp.f: 178
                  RES = ONE / Y                                       ! gamma.irp.f: 179
               ELSE                                                   ! gamma.irp.f: 180
                  RES = XINF                                          ! gamma.irp.f: 181
                  GO TO 900                                           ! gamma.irp.f: 182
            END IF                                                    ! gamma.irp.f: 183
         ELSE IF (Y .LT. TWELVE) THEN                                 ! gamma.irp.f: 184
            Y1 = Y                                                    ! gamma.irp.f: 185
            IF (Y .LT. ONE) THEN                                      ! gamma.irp.f: 186
                  Z = Y                                               ! gamma.irp.f: 190
                  Y = Y + ONE                                         ! gamma.irp.f: 191
               ELSE                                                   ! gamma.irp.f: 192
                  N = INT(Y) - 1                                      ! gamma.irp.f: 196
                  Y = Y -      N                                      ! gamma.irp.f: 197
                  Z = Y - ONE                                         ! gamma.irp.f: 198
            END IF                                                    ! gamma.irp.f: 199
            XNUM = ZERO                                               ! gamma.irp.f: 203
            XDEN = ONE                                                ! gamma.irp.f: 204
DO I = 1, 8                                                           ! gamma.irp.f: 205
               XNUM = (XNUM + P(I)) * Z                               ! gamma.irp.f: 206
               XDEN = XDEN * Z + Q(I)                                 ! gamma.irp.f: 207
  enddo                                                               ! gamma.irp.f: 208
            RES = XNUM / XDEN + ONE                                   ! gamma.irp.f: 209
            IF (Y1 .LT. Y) THEN                                       ! gamma.irp.f: 210
                  RES = RES / Y1                                      ! gamma.irp.f: 214
               ELSE IF (Y1 .GT. Y) THEN                               ! gamma.irp.f: 215
DO I = 1, N                                                           ! gamma.irp.f: 219
                     RES = RES * Y                                    ! gamma.irp.f: 220
                     Y = Y + ONE                                      ! gamma.irp.f: 221
  enddo                                                               ! gamma.irp.f: 222
            END IF                                                    ! gamma.irp.f: 223
         ELSE                                                         ! gamma.irp.f: 224
            IF (Y .LE. XBIG) THEN                                     ! gamma.irp.f: 228
                  YSQ = Y * Y                                         ! gamma.irp.f: 229
                  SUM = C(7)                                          ! gamma.irp.f: 230
DO I = 1, 6                                                           ! gamma.irp.f: 231
                     SUM = SUM / YSQ + C(I)                           ! gamma.irp.f: 232
  enddo                                                               ! gamma.irp.f: 233
                  SUM = SUM/Y - Y + SQRTPI                            ! gamma.irp.f: 234
                  SUM = SUM + (Y-HALF)*LOG(Y)                         ! gamma.irp.f: 235
                  RES = EXP(SUM)                                      ! gamma.irp.f: 236
               ELSE                                                   ! gamma.irp.f: 237
                  RES = XINF                                          ! gamma.irp.f: 238
                  GO TO 900                                           ! gamma.irp.f: 239
            END IF                                                    ! gamma.irp.f: 240
      END IF                                                          ! gamma.irp.f: 241
      IF (PARITY) then                                                ! gamma.irp.f: 245
 RES = -RES                                                           ! gamma.irp.f: 245
  endif                                                               ! gamma.irp.f: 245
      IF (FACT .NE. ONE) then                                         ! gamma.irp.f: 246
 RES = FACT / RES                                                     ! gamma.irp.f: 246
  endif                                                               ! gamma.irp.f: 246
  900 DGAMMA = RES                                                    ! gamma.irp.f: 248
      RETURN                                                          ! gamma.irp.f: 249
end                                                                   ! gamma.irp.f: 251
