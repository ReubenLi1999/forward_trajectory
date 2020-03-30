subroutine interplotation(UTCMJD,xp,yp,dx,dy,DUT1)
    !����
    !UTCMJD   UTCʱ��ļ�Լ�����գ������ڼ����ļ���ȷ��λ������ֵ
    !���
    !xp,yp    ���Ƶ�������������ֵ����������λΪ����
    !dx,dy    CIP���������������ֵ����������λΪ����
    implicit none
    !������֪���Ƴ���
    integer(kind=8)ss
    integer(kind=8)::j, flag = 0
    double precision num
    character(50) str
    !����������
    double precision UTCMJD
    !���������
    double precision xp,yp,dx,dy,DUT1!xp,yp,dx,dy,DUT1�ǲ�ֵ����
    !������֪��
    REAL(KIND=16),ALLOCATABLE::MJD(:),x(:),y(:),xx(:),yy(:),UT1UTC(:)!x,y��ԭʼ���ƣ�x,y��IERS��xx,yy��ԭʼ��dx,dy��IERS
    double precision DAS2R
    PARAMETER (DAS2R=4.848136811095359935899141D-6)!�����뻯Ϊ����
    !����Ŀ���ļ�
    open(11,file='input/jiyi.dat')
    ss=0
    do
        read(11,*, iostat=flag)
        If( flag /= 0 ) Exit
        ss=ss+1
    end do
    close(11)
    ALLOCATE(MJD(SS),X(SS),Y(SS),xx(ss),yy(ss),UT1UTC(ss))
    open(111,file='input/jiyi.dat') !���ļ���ʽΪ DATE     MJD       x       y      UT1-UTC      dX     dY
                              !           (0 h UTC)            mas     mas       ms         mas    mas ���޵�ͷ��
    do j=1,ss
        read(111,*) num,num,num,MJD(j),x(j),y(j),UT1UTC(J),xx(j),yy(j)
        x(j)=x(j)/1000d0
        y(j)=y(j)/1000d0
        UT1UTC(J)=UT1UTC(J)/1000D0
        xx(j)=xx(j)/1000d0
        yy(j)=yy(j)/1000d0
    end do
    do j=2,ss-1

        if(UTCMJD>MJD(j-1).AND.UTCMJD<MJD(j+1))then
            xp=x(j-1)*((UTCMJD-MJD(j))*(UTCMJD-MJD(j+1)))/((MJD(j-1)-MJD(j))*(MJD(j-1)-MJD(j+1)))&
            +x(j)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j+1)))/((MJD(j)-MJD(j-1))*(MJD(j)-MJD(j+1)))&
            +x(j+1)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j)))/((MJD(j+1)-MJD(j-1))*(MJD(j+1)-MJD(j)))
            yp=y(j-1)*((UTCMJD-MJD(j))*(UTCMJD-MJD(j+1)))/((MJD(j-1)-MJD(j))*(MJD(j-1)-MJD(j+1)))&
            +y(j)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j+1)))/((MJD(j)-MJD(j-1))*(MJD(j)-MJD(j+1)))&
            +y(j+1)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j)))/((MJD(j+1)-MJD(j-1))*(MJD(j+1)-MJD(j)))
            xp=xp*DAS2R
            yp=yp*DAS2R

            DUT1=UT1UTC(j-1)*((UTCMJD-MJD(j))*(UTCMJD-MJD(j+1)))/((MJD(j-1)-MJD(j))*(MJD(j-1)-MJD(j+1)))&
            +UT1UTC(j)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j+1)))/((MJD(j)-MJD(j-1))*(MJD(j)-MJD(j+1)))&
            +UT1UTC(j+1)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j)))/((MJD(j+1)-MJD(j-1))*(MJD(j+1)-MJD(j)))

            dx=xx(j-1)*((UTCMJD-MJD(j))*(UTCMJD-MJD(j+1)))/((MJD(j-1)-MJD(j))*(MJD(j-1)-MJD(j+1)))&
            +xx(j)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j+1)))/((MJD(j)-MJD(j-1))*(MJD(j)-MJD(j+1)))&
            +xx(j+1)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j)))/((MJD(j+1)-MJD(j-1))*(MJD(j+1)-MJD(j)))
            dy=yy(j-1)*((UTCMJD-MJD(j))*(UTCMJD-MJD(j+1)))/((MJD(j-1)-MJD(j))*(MJD(j-1)-MJD(j+1)))&
            +yy(j)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j+1)))/((MJD(j)-MJD(j-1))*(MJD(j)-MJD(j+1)))&
            +yy(j+1)*((UTCMJD-MJD(j-1))*(UTCMJD-MJD(j)))/((MJD(j+1)-MJD(j-1))*(MJD(j+1)-MJD(j)))
            dx=dx*DAS2R
            dy=dy*DAS2R

        end if
    end do
    close(111)
    DEALLOCATE(MJD,X,Y,xx,yy,UT1UTC)
    end subroutine
