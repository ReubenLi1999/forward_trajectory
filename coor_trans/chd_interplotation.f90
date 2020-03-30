subroutine interplotation(UTCMJD,xp,yp,dx,dy,DUT1)
    !输入
    !UTCMJD   UTC时间的简约儒略日，方便在极移文件中确定位置做插值
    !输出
    !xp,yp    极移的两个分量，插值结果，输出单位为弧度
    !dx,dy    CIP坐标的修正量，插值结果，输出单位为弧度
    implicit none
    !定义已知极移长度
    integer(kind=8)ss
    integer(kind=8)::j, flag = 0
    double precision num
    character(50) str
    !定义输入量
    double precision UTCMJD
    !定义输出量
    double precision xp,yp,dx,dy,DUT1!xp,yp,dx,dy,DUT1是插值后结果
    !定义已知量
    REAL(KIND=16),ALLOCATABLE::MJD(:),x(:),y(:),xx(:),yy(:),UT1UTC(:)!x,y是原始极移（x,y）IERS，xx,yy是原始（dx,dy）IERS
    double precision DAS2R
    PARAMETER (DAS2R=4.848136811095359935899141D-6)!将角秒化为弧度
    !定义目标文件
    open(11,file='input/jiyi.dat')
    ss=0
    do
        read(11,*, iostat=flag)
        If( flag /= 0 ) Exit
        ss=ss+1
    end do
    close(11)
    ALLOCATE(MJD(SS),X(SS),Y(SS),xx(ss),yy(ss),UT1UTC(ss))
    open(111,file='input/jiyi.dat') !该文件格式为 DATE     MJD       x       y      UT1-UTC      dX     dY
                              !           (0 h UTC)            mas     mas       ms         mas    mas （无道头）
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
