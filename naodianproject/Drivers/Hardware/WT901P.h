#ifndef __WT901P_H
#define __WT901P_H

void WT901_Init(void);
void RSW_control(void);
void change_baudrate(void);
void sendSensorData(float ax, float ay, float az, 
                    float wx, float wy, float wz,
                    float Roll, float Pitch, float Yaw,
					float w,float,float,float) ;
void Data_processing(uint8_t* RxData,uint16_t Size);
#define g 9.8f
#define DT 0.1f
void zero_offset(float *a_bias,float world_ax,float world_ay,float world_az);
void zero_offset(float *a_bias,float world_ax,float world_ay,float world_az);
void RRATE_control(void);
#endif
