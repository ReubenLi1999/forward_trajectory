subroutine ITRSandGCRS_P(flag,DATE1,DATE2,LEAPSECOND,P1,P2)
    !-----------本程序支持2020年1月的数据--------------------!
    implicit none
    !--------------转换方向------------------
    integer flag                                      !flag=1，表示GCRS向ITRS转换；flag=2，表示ITRS向GCRS转换
    
    !--------------时间----------------------
    DOUBLE PRECISION date1,date2  !date1是儒略日（TT）天，xxx天，date2是儒略日（TT）超出整天部分，xxx天：e.g.2020年1月1日中午12h(TT)，date1=2458849.5d0，date2=0.5d0
    DOUBLE PRECISION DJ1,DJ2                         !UT1儒略日的两部分
    
    !--------------坐标和速度----------------
    DOUBLE PRECISION Px1,Py1,Pz1 !输入的位置
    DOUBLE PRECISION P1(3,1)
    DOUBLE PRECISION Px2,Py2,Pz2 !输出的位置
    DOUBLE PRECISION P2(3,1)
    
    !--------------转换矩阵------------------
    DOUBLE PRECISION RPOM(3,3),RBPN(3,3),RT2C(3,3)  !极移矩阵，岁差章动矩阵，ITRS到GCRS的转换矩阵
    DOUBLE PRECISION RC2T(3,3) !GCRS到ITRS的转换矩阵
    DOUBLE PRECISION RC2I(3,3) !GCRS到CIRS的转换矩阵,GCRS到TIRS的转换矩阵
    DOUBLE PRECISION RB(3,3), RP(3,3), RBP(3,3),RN(3,3)!
    
    !--------------转换所需参数--------------
    DOUBLE PRECISION iau_GST06,iau_SP00!计算GAST\sp的函数
    DOUBLE PRECISION DPSI, DEPS, EPSA!章动角，章动角修正值
    DOUBLE PRECISION GST!格林尼治真恒星时
    DOUBLE PRECISION sp                 !极移矩阵中的一个量
    DOUBLE PRECISION X, Y, S            !CIP坐标，“CIP定位器”（一个数值）
    
    !--------------数据本身部分信息----------
    INTEGER Leapsecond   !可由IERS公告得到
    DOUBLE PRECISION UTCMJD,xp,yp,dx,dy !UTC简约儒略日，极移的两个坐标，极移的修正量
    DOUBLE PRECISION DUT1!UT1-UTC=DUT1
    
    !---------------转换开始-------------------
    
    UTCMJD=DATE1+DATE2-(32.184D0+LEAPSECOND)/86400D0-2400000.5D0
    call interplotation(UTCMJD,xp,yp,dx,dy,DUT1)!!!!!!!!!!这是自己的程序
    
    call iau_NUT06A ( DATE1, DATE2, DPSI, DEPS )
    call iau_PN06 ( DATE1, DATE2, DPSI, DEPS,EPSA, RB, RP, RBP, RN, RBPN )
        
    DJ1=DATE1
    DJ2=DATE2-(LEAPSECOND+32.184d0-DUT1)/86400D0
    GST=iau_GST06 ( DJ1, DJ2, date1, date2, RBPN )
    call iau_RZ ( GST, RBPN )
        
    sp=iau_SP00 ( DATE1, DATE2 )
    call iau_POM00 ( XP, YP, SP, RPOM )
        
    call iau_RXR ( RPOM, RBPN, RC2T )
        

    IF(flag==1)THEN                               !!!!flag=1，表示GCRS向ITRS转换
        P2=MATMUL(RC2T,P1)
    ELSE IF(flag==2)THEN                          !!!!flag=2，表示ITRS向GCRS转换
        call iau_TR (RC2T,RT2C)
        P2=MATMUL(RT2C,P1)
    END IF
    
end subroutine ITRSandGCRS_P