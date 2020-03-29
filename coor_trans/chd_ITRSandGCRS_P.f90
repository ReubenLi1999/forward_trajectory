subroutine ITRSandGCRS_P(flag,DATE1,DATE2,LEAPSECOND,P1,P2)
    !-----------������֧��2020��1�µ�����--------------------!
    implicit none
    !--------------ת������------------------
    integer flag                                      !flag=1����ʾGCRS��ITRSת����flag=2����ʾITRS��GCRSת��
    
    !--------------ʱ��----------------------
    DOUBLE PRECISION date1,date2  !date1�������գ�TT���죬xxx�죬date2�������գ�TT���������첿�֣�xxx�죺e.g.2020��1��1������12h(TT)��date1=2458849.5d0��date2=0.5d0
    DOUBLE PRECISION DJ1,DJ2                         !UT1�����յ�������
    
    !--------------������ٶ�----------------
    DOUBLE PRECISION Px1,Py1,Pz1 !�����λ��
    DOUBLE PRECISION P1(3,1)
    DOUBLE PRECISION Px2,Py2,Pz2 !�����λ��
    DOUBLE PRECISION P2(3,1)
    
    !--------------ת������------------------
    DOUBLE PRECISION RPOM(3,3),RBPN(3,3),RT2C(3,3)  !���ƾ�������¶�����ITRS��GCRS��ת������
    DOUBLE PRECISION RC2T(3,3) !GCRS��ITRS��ת������
    DOUBLE PRECISION RC2I(3,3) !GCRS��CIRS��ת������,GCRS��TIRS��ת������
    DOUBLE PRECISION RB(3,3), RP(3,3), RBP(3,3),RN(3,3)!
    
    !--------------ת���������--------------
    DOUBLE PRECISION iau_GST06,iau_SP00!����GAST\sp�ĺ���
    DOUBLE PRECISION DPSI, DEPS, EPSA!�¶��ǣ��¶�������ֵ
    DOUBLE PRECISION GST!�������������ʱ
    DOUBLE PRECISION sp                 !���ƾ����е�һ����
    DOUBLE PRECISION X, Y, S            !CIP���꣬��CIP��λ������һ����ֵ��
    
    !--------------���ݱ�������Ϣ----------
    INTEGER Leapsecond   !����IERS����õ�
    DOUBLE PRECISION UTCMJD,xp,yp,dx,dy !UTC��Լ�����գ����Ƶ��������꣬���Ƶ�������
    DOUBLE PRECISION DUT1!UT1-UTC=DUT1
    
    !---------------ת����ʼ-------------------
    
    UTCMJD=DATE1+DATE2-(32.184D0+LEAPSECOND)/86400D0-2400000.5D0
    call interplotation(UTCMJD,xp,yp,dx,dy,DUT1)!!!!!!!!!!�����Լ��ĳ���
    
    call iau_NUT06A ( DATE1, DATE2, DPSI, DEPS )
    call iau_PN06 ( DATE1, DATE2, DPSI, DEPS,EPSA, RB, RP, RBP, RN, RBPN )
        
    DJ1=DATE1
    DJ2=DATE2-(LEAPSECOND+32.184d0-DUT1)/86400D0
    GST=iau_GST06 ( DJ1, DJ2, date1, date2, RBPN )
    call iau_RZ ( GST, RBPN )
        
    sp=iau_SP00 ( DATE1, DATE2 )
    call iau_POM00 ( XP, YP, SP, RPOM )
        
    call iau_RXR ( RPOM, RBPN, RC2T )
        

    IF(flag==1)THEN                               !!!!flag=1����ʾGCRS��ITRSת��
        P2=MATMUL(RC2T,P1)
    ELSE IF(flag==2)THEN                          !!!!flag=2����ʾITRS��GCRSת��
        call iau_TR (RC2T,RT2C)
        P2=MATMUL(RT2C,P1)
    END IF
    
end subroutine ITRSandGCRS_P