#include "pwm_uart.hpp"
#include "ekf_sensor.hpp"


//e エレベータ a エルロン r ラダー t スロットル
float h=0.01;
float kpe,kie,kde;
float kpa,kia,kda;
float kpr,kir,kdr;
float shigma_e,shiguma_a,shiguma_r;
float ref_e,ref_r,ref_a,err_e,err_a,err_r;
float sk_e,sk_a,sk_r;
float olderr_e,olderr_a,olderr_r;
float dk_a,dk_e,dk_r;
float se,sa,sr;
/* Private macro -------------------------------------------------------------*/

int main(void)
{
  Matrix<float, 7 ,1> xp = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> xe = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> x_sim = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
  Matrix<float, 6 ,1> z = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_sim = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_noise = MatrixXf::Zero(6,1);
  Matrix<float, 3, 1> omega_m = MatrixXf::Zero(3, 1);
  Matrix<float, 3, 1> omega_sim;
  Matrix<float, 3, 1> domega;
  Matrix<float, 3, 1> domega_sim;
  Matrix<float, 3, 3> Q = MatrixXf::Identity(3, 3)*1;
  Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
  Matrix<float, 7 ,3> G;
  Matrix<float, 3 ,1> beta;
  float t=0.0,dt=0.01;
  uint64_t s_time=0,e_time=0,d_time=0;
  double phl,theta,psi;
  short i,waittime=5;
  float p_com[10]={0.1*PI, 0, 0.01*PI, 0, -0.01*PI, 0, 0.02*PI, 0, -0.05*PI, 0};
  float q_com[10]={0.1*PI, 0, 0      , 0, -0.01*PI, 0, 0.02*PI, 0, -0.02*PI, 0};
  float r_com[10]={0.1*PI, 0,-0.02*PI, 0,      -PI, 0, 0.02*PI, 0,  0.02*PI, 0};
  float endtime=10000.0;
  float control_period=5.0;
  float control_time=5.0;
  int counter=0;
  int sample=1;
  int control_counter=0;
  int control_counter_max=0;
  double pi=3.14159265358;
  float ax,ay,az,wp,wq,wr,mx,my,mz; 
  float p11=-0.60526553,p12=0.79021444,p13=-0.09599364,p21=0.78892428,p22=0.61155945,p23=0.05994594,p31=-0.10607597,p32=0.0394485,p33=0.99357521;
  float r1=-4.96008102e-6,r2=-4.60371715e-6,r3=-4.17649591e-6;
  float w;
  float dmx,dmy,dmz,dxx,dxy,dxz;
  float ddxx, ddxy,ddxz;
  float dddxx,dddxy,dddxz;
  float ddddxx,ddddxy,ddddxz;
 //pwm_uart
  float duty_rr,duty_rl,duty_fl,duty_fr; 
  float Duty_rr,Duty_rl,Duty_fl,Duty_fr;
  float seigyo=0.5;  

  
  xe << 1.0, 0.0, 0.0, 0.0, -0.078, 0.0016, 0.00063;
  xp =xe;

  G <<  0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        1.0,0.0,0.0, 
        0.0,1.0,0.0, 
        0.0,0.0,1.0;

  beta << 0.003, 0.003, 0.003;

  P <<  1,0,0,0,0,0,0,  
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,  
        0,0,0,1,0,0,0, 
        0,0,0,0,1,0,0,  
        0,0,0,0,0,1,0,  
        0,0,0,0,0,0,1;

  Q << 7.34944e-6,0,0,
       0,6.861e-6,0,
       0,0,5.195e-6;

  R << 3.608e-6,0,0,0,0,0,
       0,6.261e-6,0,0,0,0,
       0,0,1.889e-5,0,0,0,
       0,0,0,1.0,0,0,
       0,0,0,0,1.0,0,
       0,0,0,0,0,1.0;


  stdio_init_all();
  pwm_settei();
  serial_settei();
 // imu_mag_init();
  
  for (i=0;i<waittime;i++)
  {
    printf("#Please wait %d[s] ! \n",waittime-i);
    sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
  
  while(1/*t<endtime*/)
  {
   /* s_time=time_us_64();
    sleep_ms(10);
    //Control
    if(t>control_time)
    {
      control_time = control_time + control_period;
      control_counter++;
      if(control_counter>control_counter_max)control_counter=0;
    }
    ax=   -acceleration_mg[0]*1/1000*GRAV;
    ay=   -acceleration_mg[1]*1/1000*GRAV;
    az=    acceleration_mg[2]*1/1000*GRAV;
    wp=    angular_rate_mdps[0]/1000*pi/180;
    wq=    angular_rate_mdps[1]/1000*pi/180;
    wr=   -angular_rate_mdps[2]/1000*pi/180;
    dmx=  -(magnetic_field_mgauss[0]-310);
    dmy=   (magnetic_field_mgauss[1]-10);
    dmz=  -(magnetic_field_mgauss[2]);
 // elevator();    
    //校正作業
    dxx=p11*dmx+p21*dmy+p31*dmz;
    dxy=p12*dmx+p22*dmy+p32*dmz;
    dxz=p13*dmx+p23*dmy+p33*dmz;

    ddxx=dxx-22.760831749415342;
    ddxy=dxy-19.734355196006327;
    ddxz=dxz-141.33565570453044;


    w=-1.087745370038146;
  
    dddxx=ddxx*0.0020572671658147883;
    dddxy=ddxy*0.0021354074993493823;
    dddxz=ddxz*0.0019594870993397107;

    mx=p11*dddxx+p22*dddxy+p13*dddxz;
    my=p21*dddxx+p22*dddxy+p23*dddxz;
    mz=p31*dddxx+p32*dddxy+p33*dddxz;
 
    omega_m <<wp,wq,wr;
    z       <<ax,ay,az,mx,my,mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    ekf(xp, xe, P, z, omega_m, Q, R, G*dt, beta, dt);
    if(counter%sample==0)
    {
		imu_mag_data_read();
                phl=Phl(xe);
                theta=Theta(xe);
                psi=Psi(xe);
        
    
   }*/ 
/*    //姿勢安定化
    //最大角度30°
    ref_e=Data2 * 0.523598775;
    ref_a=Data4 * 0.523598775;
    ref_r=Data1 * 0.523598775;
    //エレベータ
    olderr_e=err_e;
    err_e=ref_e - psi;
    if (sk_e<=30000){
      sk_e=sk_e+olderr_e;
    }
    else if(-30000<=sk_e){
      sk_e=sk_e+olderr_e;
    }
    dk_e=(err_e-olderr_e)*100;
    se=kpe*err_e+kie*sk_e+kde*dk_e;

    //エルロン
 
    olderr_a=err_a;
    err_a=ref_a - psi;
    if (sk_a<=30000){
      sk_a=sk_a+olderr_a;
    }
    else if(-30000<=sk_a){
      sk_a=sk_a+olderr_a;
    }
    dk_a=(err_a-olderr_a)*100;
    sa=kpa*err_a+kia*sk_a+kda*dk_a;

    //ラダー
        
    olderr_r=err_r;
    err_r=ref_r - psi;
    if (sk_r<=30000){
      sk_r=sk_r+olderr_r;
    }
    else if(-30000<=sk_r){
      sk_r=sk_r+olderr_r;
    }
    dk_r=(err_r-olderr_r)*100;
    sr=kpr*err_r+kir*sk_r+kdr*dk_r;


    Duty_fr=Data3+se-sa-sr;
    Duty_fl=Data3+se+sa+sr;
    Duty_rr=Data3-se-sa+sr;         
    Duty_rl=Data3-se+sa-sr;
*/
   
    Duty_fr=Data3;
    Duty_fl=Data3;
    Duty_rr=Data3;         
    Duty_rl=Data3;

    tight_loop_contents();
    duty_rr=(float)(DUTYMAX-DUTYMIN)*Duty_rr+DUTYMIN;
    duty_fr=(float)(DUTYMAX-DUTYMIN)*Duty_fr+DUTYMIN;
    duty_rl=(float)(DUTYMAX-DUTYMIN)*Duty_rl+DUTYMIN;
    duty_fl=(float)(DUTYMAX-DUTYMIN)*Duty_fl+DUTYMIN;
    if (duty_rr>DUTYMAX)duty_rr=DUTYMAX;
    if (duty_rr<DUTYMIN)duty_rr=DUTYMIN;
    if (duty_fr>DUTYMAX)duty_fr=DUTYMAX;
    if (duty_fr>DUTYMAX)duty_fr=DUTYMAX;
    if (duty_rl>DUTYMAX)duty_rl=DUTYMAX;
    if (duty_rl<DUTYMIN)duty_rl=DUTYMIN;
    if (duty_fl<DUTYMIN)duty_fl=DUTYMIN;
    if (duty_fl<DUTYMIN)duty_fl=DUTYMIN;
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, duty_rr);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, duty_fr);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, duty_rl);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, duty_fl);
    //printf("%04f %04f %04f %04f %04f %04f %04f %04f %04f %04f \n",Data1,Data2,Data3,Data4,Data5,Data6,xe(0,0), xe(1,0), xe(2,0),xe(3,0));  
    sleep_ms(10);
  }
}



#if 0
int main(void)
{
  Matrix<float, 7 ,1> xp = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> xe = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,1> x_sim = MatrixXf::Zero(7,1);
  Matrix<float, 7 ,7> P = MatrixXf::Identity(7,7);
  Matrix<float, 6 ,1> z = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_sim = MatrixXf::Zero(6,1);
  Matrix<float, 6 ,1> z_noise = MatrixXf::Zero(6,1);
  Matrix<float, 3, 1> omega_m = MatrixXf::Zero(3, 1);
  Matrix<float, 3, 1> omega_sim;
  Matrix<float, 3, 1> domega;
  Matrix<float, 3, 1> domega_sim;
  Matrix<float, 3, 3> Q = MatrixXf::Identity(3, 3)*1;
  Matrix<float, 6, 6> R = MatrixXf::Identity(6, 6)*1;
  Matrix<float, 7 ,3> G;
  Matrix<float, 3 ,1> beta;
  float t=0.0,dt=0.01;
  uint64_t s_time=0,e_time=0,d_time=0;
  double phl,theta,psi;
  short i,waittime=5;
  float p_com[10]={0.1*PI, 0, 0.01*PI, 0, -0.01*PI, 0, 0.02*PI, 0, -0.05*PI, 0};
  float q_com[10]={0.1*PI, 0, 0      , 0, -0.01*PI, 0, 0.02*PI, 0, -0.02*PI, 0};
  float r_com[10]={0.1*PI, 0,-0.02*PI, 0,      -PI, 0, 0.02*PI, 0,  0.02*PI, 0};
  float endtime=10000.0;
  float control_period=5.0;
  float control_time=5.0;
  int counter=0;
  int sample=1;
  int control_counter=0;
  int control_counter_max=0;
  double pi=3.14159265358;
  float ax,ay,az,wp,wq,wr,mx,my,mz; 
  float p11=-0.60526553,p12=0.79021444,p13=-0.09599364,p21=0.78892428,p22=0.61155945,p23=0.05994594,p31=-0.10607597,p32=0.0394485,p33=0.99357521;
  float r1=-4.96008102e-6,r2=-4.60371715e-6,r3=-4.17649591e-6;
  float w;
  float dmx,dmy,dmz,dxx,dxy,dxz;
  float ddxx, ddxy,ddxz;
  float dddxx,dddxy,dddxz;
  float ddddxx,ddddxy,ddddxz;
 //pwm_uart
  float duty_rr,duty_rl,duty_fl,duty_fr;
  float seigyo=0.5;  
  stdio_init_all();
  pwm_settei();
  serial_settei();
  
 // elevator();
  sleep_ms(2000);
  //Variable Initalize
  xe << 1.0, 0.0, 0.0, 0.0, -0.078, 0.0016, 0.00063;
  xp =xe;
  //x_sim << 1.0, 0.0, 0.0, 0.0, 0.01, 0.02, 0.03;
 // observation_equation(x_sim, z_sim, GRAV, MN, MD);

  G <<  0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        0.0,0.0,0.0, 
        1.0,0.0,0.0, 
        0.0,1.0,0.0, 
        0.0,0.0,1.0;

  beta << 0.003, 0.003, 0.003;

  P <<  1,0,0,0,0,0,0,  
        0,1,0,0,0,0,0,
        0,0,1,0,0,0,0,  
        0,0,0,1,0,0,0, 
        0,0,0,0,1,0,0,  
        0,0,0,0,0,1,0,  
        0,0,0,0,0,0,1;

  Q << 7.34944e-6,0,0,
       0,6.861e-6,0,
       0,0,5.195e-6;

  R << 3.608e-6,0,0,0,0,0,
       0,6.261e-6,0,0,0,0,
       0,0,1.889e-5,0,0,0,
       0,0,0,1.0,0,0,
       0,0,0,0,1.0,0,
       0,0,0,0,0,1.0;
  
  //Initilize Console Input&Output
  imu_mag_init();

  


  //Start up wait for Pico
  for (i=0;i<waittime;i++)
  {
    printf("#Please wait %d[s] ! \n",waittime-i);
    sleep_ms(1000);
  }
  printf("#Start Kalman Filter\n");
    
  while(t<endtime)
  {
    s_time=time_us_64();
    sleep_ms(10);
    //Control
    if(t>control_time)
    {
      control_time = control_time + control_period;
      control_counter++;
      if(control_counter>control_counter_max)control_counter=0;
    }
    ax=   -acceleration_mg[0]*1/1000*GRAV;
    ay=   -acceleration_mg[1]*1/1000*GRAV;
    az=    acceleration_mg[2]*1/1000*GRAV;
    wp=    angular_rate_mdps[0]/1000*pi/180;
    wq=    angular_rate_mdps[1]/1000*pi/180;
    wr=   -angular_rate_mdps[2]/1000*pi/180;
    dmx=  -(magnetic_field_mgauss[0]-310);
    dmy=   (magnetic_field_mgauss[1]-10);
    dmz=  -(magnetic_field_mgauss[2]);




    //校正作業
    dxx=p11*dmx+p21*dmy+p31*dmz;
    dxy=p12*dmx+p22*dmy+p32*dmz;
    dxz=p13*dmx+p23*dmy+p33*dmz;

    ddxx=dxx-22.760831749415342;
    ddxy=dxy-19.734355196006327;
    ddxz=dxz-141.33565570453044;


    w=-1.087745370038146;
  
    dddxx=ddxx*0.0020572671658147883;
    dddxy=ddxy*0.0021354074993493823;
    dddxz=ddxz*0.0019594870993397107;

    mx=p11*dddxx+p22*dddxy+p13*dddxz;
    my=p21*dddxx+p22*dddxy+p23*dddxz;
    mz=p31*dddxx+p32*dddxy+p33*dddxz;
 
    omega_m <<wp,wq,wr;
    z       <<ax,ay,az,mx,my,mz;//ここに入れる
    //--Begin Extended Kalman Filter--
    ekf(xp, xe, P, z, omega_m, Q, R, G*dt, beta, dt);
   // e_time=time_us_64();
    //--End   Extended Kalman Filter--
    printf("%9.2f %9.6f %9.6f %9.6f \n",t,mx,my,mz);

    Duty_fr=Data3+((Data2-.5)-(Data4-.5)-(Data1-.5))*seigyo;
    Duty_fl=Data3+(-1.5+Data2+Data4+Data1)*seigyo;
    Duty_rr=Data3+(-(Data2-.5)-(Data4-.5)+(Data1-.5))*seigyo;         
    Duty_rl=Data3+(-(Data2-.5)+(Data4-.5)-(Data1-.5))*seigyo;

    tight_loop_contents();
    duty_rr=(float)(DUTYMAX-DUTYMIN)*Duty_rr+DUTYMIN;
    duty_fr=(float)(DUTYMAX-DUTYMIN)*Duty_fr+DUTYMIN;
    duty_rl=(float)(DUTYMAX-DUTYMIN)*Duty_rl+DUTYMIN;
    duty_fl=(float)(DUTYMAX-DUTYMIN)*Duty_fl+DUTYMIN;
    if (duty_rr>DUTYMAX)duty_rr=DUTYMAX;
    if (duty_rr<DUTYMIN)duty_rr=DUTYMIN;
    if (duty_fr>DUTYMAX)duty_fr=DUTYMAX;
    if (duty_fr>DUTYMAX)duty_fr=DUTYMAX;
    if (duty_rl>DUTYMAX)duty_rl=DUTYMAX;
    if (duty_rl<DUTYMIN)duty_rl=DUTYMIN;
    if (duty_fl<DUTYMIN)duty_fl=DUTYMIN;
    if (duty_fl<DUTYMIN)duty_fl=DUTYMIN;
    pwm_set_chan_level(slice_num[0], PWM_CHAN_A, duty_rr);
    pwm_set_chan_level(slice_num[0], PWM_CHAN_B, duty_fr);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_A, duty_rl);
    pwm_set_chan_level(slice_num[1], PWM_CHAN_B, duty_fl);
    printf("%04f %04f %04f %04f %04f %04f %04f %04f %04f %04f \n",Data1,Data2,Data3,Data4,Data5,Data6,duty_rr,duty_rl,duty_fl,duty_fr);  
     sleep_ms(10);
    //Result output
    if(counter%sample==0)
    {
		imu_mag_data_read();
                phl=Phl(xe);
                theta=Theta(xe);
                psi=Psi(xe);
              //  printf("%9.6f %9.6f %9.6f %9.6f  ",t,phl,theta,psi);
     /* printf("%9.2f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n", 
               t, 	xe(0,0), xe(1,0), xe(2,0),xe(3,0), xe(4,0), xe(5,0),xe(6,0)
     // printf("%9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f %9.6f\n",ax,ay,az,wp,wq,wr,mx,my,mz);       
	 x_sim(0,0), x_sim(1,0), x_sim(2,0), x_sim(3,0),
	        x_sim(4,0), x_sim(5,0), x_sim(6,0),
                p_com[control_counter], q_com[control_counter], r_com[control_counter],
                e_time-s_time);  
      printf( "%9.6fIMU-[mg]:\t%9.6f\t%9.6f\t%9.6f\t[mdps]:\t%9.6f\t%9.6f\t%9.6f\t",t
              ,acceleration_mg[0]*1/1000*GRAV, acceleration_mg[1]*1/1000*GRAV, acceleration_mg[2]*1/1000*GRAV*(-1),
                angular_rate_mdps[0]/1000*pi/180, angular_rate_mdps[1]/1000*pi/180, angular_rate_mdps[2]/1000*pi/180*(-1));
      printf( "MAG-[mG]:\t%9.6f\t%9.6f\t%9.6f\r\n"
               ,(magnetic_field_mgauss[0]-310)*(-1), magnetic_field_mgauss[1]-10,magnetic_field_mgauss[2]*(-1));*/
                //Control
               // domega<<xe(4,0), xe(5,0), xe(6,0);
              	//printf("%f,%f,%f\n",angular_rate_mdps[0], angular_rate_mdps[1], angular_rate_mdps[2]*(-1));*/
                //omega=omega_sim+domega;
                //Simulation
               // observation_equation(xe, z, GRAV, MN, MD);
                //z=z_sim;
                //rk4(quatdot, t, dt, quat_sim, omega_sim);
                //x_sim << quat_sim(0,0), quat_sim(1,0), quat_sim(2,0), quat_sim(3,0), 0,0,0;
                //t=t+dt;
                //Begin Extended Kalman Filter
                // s_time=time_us_64();
                //ekf(xe, P, z, omega_m, beta, dt);
      printf("%9.2f %9.6f %9.6f %9.6f \n",t,mx,my,mz);

     }
 
    counter++;
    
      	
    
    t=t+dt;
    while (time_us_64()-s_time<10000);

  }
   return 0;

}
#endif
