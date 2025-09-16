#include "stm32f4xx.h"                  // Device header
#include "usart.h"
#include "stdio.h"
#include "WT901P.h"
#include "matrix.h"
#include "math.h"
void WT901_accInit(void)
{
	printf("Please place it horizontally\r\n");
	for(uint16_t i=1;i<=10;i++)
	{
	printf("%d\r\n",i);
	HAL_Delay(1000);}
	uint8_t txData1[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData2[5]={0xFF,0xAA,0x01,0x01,0x00};
	HAL_UART_Transmit(&huart2,txData2,5,HAL_MAX_DELAY);
	HAL_Delay(4000);//加速度校准
	
	uint8_t txData3[5]={0xFF,0xAA,0x01,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData3,5,HAL_MAX_DELAY);
	HAL_Delay(100);//退出校准

	uint8_t txData4[5]={0xFF,0xAA,0x00,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData4,5,HAL_MAX_DELAY);//保存
}

void WT901_angleInit(void)
{
	uint8_t txData1[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData2[5]={0xFF,0xAA,0x01,0x08,0x00};
	HAL_UART_Transmit(&huart2,txData2,5,HAL_MAX_DELAY);
	HAL_Delay(3000);//角度校准

	uint8_t txData3[5]={0xFF,0xAA,0x00,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData3,5,HAL_MAX_DELAY);//保存
}

void change_baudrate(void)
{
	uint8_t txData1[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData2[5]={0xFF,0xAA,0x04,0x06,0x00};
	HAL_UART_Transmit(&huart2,txData2,5,HAL_MAX_DELAY);
	
	while (__HAL_UART_GET_FLAG(&huart2, UART_FLAG_TC) == RESET);//等待发送完成
	HAL_UART_Abort(&huart2);//中止串口操作
	huart2.Instance = USART2;
	huart2.Init.BaudRate = 115200;
	huart2.Init.WordLength = UART_WORDLENGTH_8B;
	huart2.Init.StopBits = UART_STOPBITS_1;
	huart2.Init.Parity = UART_PARITY_NONE;
	huart2.Init.Mode = UART_MODE_TX_RX;
	huart2.Init.HwFlowCtl = UART_HWCONTROL_NONE;
	huart2.Init.OverSampling = UART_OVERSAMPLING_16;    if (HAL_UART_Init(&huart2) != HAL_OK) {
        Error_Handler();
		}//修改波特率为115200hz

	uint8_t txData3[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData3,5,HAL_MAX_DELAY);
	HAL_Delay(200);//重新解锁
	
	uint8_t txData4[5]={0xFF,0xAA,0x00,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData4,5,HAL_MAX_DELAY);//保存
}

void WT901_gyro_Init(void)
{
	uint8_t txData1[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData2[5]={0xFF,0xAA,0x61,0x01,0x00};
	HAL_UART_Transmit(&huart2,txData2,5,HAL_MAX_DELAY);
	HAL_Delay(3000);//陀螺仪校准

	uint8_t txData3[5]={0xFF,0xAA,0x00,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData3,5,HAL_MAX_DELAY);//保存	
}

void WT901_MF_Init(void)
{
	uint8_t txData1[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData2[5]={0xFF,0xAA,0x01,0x07,0x00};
	HAL_UART_Transmit(&huart2,txData2,5,HAL_MAX_DELAY);
	printf("To begin the magnetic field calibration, slowly rotate around the three axes.\r\n");
	HAL_Delay(10);
	for(uint8_t i=1;i<21;i++)
	{printf("%d\r\n",i);
	HAL_Delay(1000);}
	//磁场校准，缓慢绕三个轴旋转1-2圈

	uint8_t txData3[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData3,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData4[5]={0xFF,0xAA,0x01,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData4,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
	
	uint8_t txData5[5]={0xFF,0xAA,0x00,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData5,5,HAL_MAX_DELAY);//保存		
}

void unlock(void)
{
	uint8_t txData1[5]={0xFF,0xAA,0x69,0x88,0xB5};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	HAL_Delay(200);//解锁
}
void save(void)
{
	uint8_t txData5[5]={0xFF,0xAA,0x00,0x00,0x00};
	HAL_UART_Transmit(&huart2,txData5,5,HAL_MAX_DELAY);//保存	
}

void RSW_control(void)
{
	unlock();
	uint8_t txData1[5]={0xFF,0xAA,0x02,0x0E,0x02};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	save();
}

void RRATE_control(void)
{
	unlock();
	uint8_t txData1[5]={0xFF,0xAA,0x03,0x09,0x00};
	HAL_UART_Transmit(&huart2,txData1,5,HAL_MAX_DELAY);
	save();
}

void WT901_Init(void)
{
	change_baudrate();
	WT901_accInit();
	WT901_angleInit();
	WT901_gyro_Init();
	WT901_MF_Init();
	RSW_control();
	printf("Initialization complete\r\n");
}

// 定义帧结构体
#pragma pack(push, 1)  // 确保紧密内存对齐
typedef struct {
    float data[13];      // 13个浮点数：ax,ay,az,wx,wy,wz,Roll,Pitch,Yaw
    uint8_t tail[4];    // 帧尾
} SensorFrame;
#pragma pack(pop)

void sendSensorData(float ax, float ay, float az, 
                    float wx, float wy, float wz,
                    float Roll, float Pitch, float Yaw,
					float w,float x,float y,float z)
{
    // 方法1：使用HAL_UART直接发送（推荐）
    SensorFrame frame;
    
    // 填充数据
    frame.data[0] = ax;
    frame.data[1] = ay;
    frame.data[2] = az;
    frame.data[3] = wx;
    frame.data[4] = wy;
    frame.data[5] = wz;
    frame.data[6] = Roll;
    frame.data[7] = Pitch;
    frame.data[8] = Yaw;
    frame.data[9] = w;
    frame.data[10] =x;
    frame.data[11] =y;
    frame.data[12] =z;
	
    
    // 设置帧尾
    frame.tail[0] = 0x00;
    frame.tail[1] = 0x00;
    frame.tail[2] = 0x80;
    frame.tail[3] = 0x7F;
    
    // 通过UART发送整个帧
    HAL_UART_Transmit(&huart1, (uint8_t*)&frame, sizeof(frame), HAL_MAX_DELAY);

    /* 
    // 方法2：使用printf发送（需要重定向，不推荐）
    uint8_t buffer[sizeof(SensorFrame)];
    memcpy(buffer, &frame, sizeof(frame));
    for(int i=0; i<sizeof(buffer); i++) {
        printf("%c", buffer[i]);  // 逐字节输出
    }
    */
}
/**
 * @brief 将本体加速度转换到世界坐标系
 * @param ax 本体坐标系x轴加速度
 * @param ay 本体坐标系y轴加速度
 * @param az 本体坐标系z轴加速度
 * @param qw 四元数w分量
 * @param qx 四元数x分量
 * @param qy 四元数y分量
 * @param qz 四元数z分量
 * @param world_ax 世界坐标系x轴加速度输出指针
 * @param world_ay 世界坐标系y轴加速度输出指针
 * @param world_az 世界坐标系z轴加速度输出指针
 */
void convert_body_to_world_acceleration(float ax, float ay, float az, 
                                        float qw, float qx, float qy, float qz,
                                        float* world_ax, float* world_ay, float* world_az) {
    float body_accel[3] = {ax, ay, az};
    float world_accel[3];
    
    // 使用四元数旋转向量
    rotate_vector_by_quaternion(body_accel, qw, qx, qy, qz, world_accel);
    
    *world_ax = world_accel[0];
    *world_ay = world_accel[1];
    *world_az = world_accel[2];
}


	uint16_t zero_offset_flat=0;
		short raw_ax=0;
	short raw_ay=0;
	short raw_az=0;
	short raw_wx=0;
	short raw_wy=0;
	short raw_wz=0;
	short raw_Roll=0;
	short raw_Pitch=0;
	short raw_Yaw=0;
	short raw_w=0;
	short raw_x=0;
	short raw_y=0;
	short raw_z=0;
	float real_ax=0;
	float real_ay=0;
	float real_az=0;
	float real_wx=0;
	float real_wy=0;
	float real_wz=0;
	float real_Roll=0;
	float real_Pitch=0;
	float real_Yaw=0;
	float real_w=1;
	float real_x=0;
	float real_y=0;
	float real_z=0;
	float world_ax=0;
	float world_ay=0;
	float world_az=0;
	uint8_t stop_cnt=0;
	float a_bias[3] = {0, 0, 0}; // 加速度计零偏
	float v_bias[3] = {0, 0, 0}; // 加速度计零偏	
	float v_prev[3] = {0, 0, 0}; // 上一时刻速度
	float s_prev[3] = {0, 0, 0}; // 上一时刻位移
//	static float filtered_ax = 0, filtered_ay = 0, filtered_az = 0; // 用于高通滤波
	const float alpha = 0.98f; // 高通滤波系数
void Data_processing(uint8_t* RxData,uint16_t Size)
{

			int index=0;
		while(index<=Size-11)
		{
			if(RxData[index]!=0x55)
			{
				index++;
				continue;
			}
			uint8_t SUMCRC=0;
			switch(RxData[index+1])
			{
				case 0x51:
					for(int i=0;i<10;i++){
					SUMCRC+=RxData[index+i];
					}
					if(SUMCRC!=RxData[index+10]){
						index+=11;
						continue;
					}
					else{
						raw_ax=(short)((short)RxData[index+3]<<8|RxData[index+2]);
						raw_ay=(short)((short)RxData[index+5]<<8|RxData[index+4]);
						raw_az=(short)((short)RxData[index+7]<<8|RxData[index+6]);
					}
					break;
				case 0x52:
					for(int i=0;i<10;i++){
					SUMCRC+=RxData[index+i];
					}
					if(SUMCRC!=RxData[index+10]){
						index+=11;
						continue;
					}
					else{
						raw_wx=(short)((short)RxData[index+3]<<8|RxData[index+2]);
						raw_wy=(short)((short)RxData[index+5]<<8|RxData[index+4]);
						raw_wz=(short)((short)RxData[index+7]<<8|RxData[index+6]);
					}					
					break;
				case 0x53:
					for(int i=0;i<10;i++){
					SUMCRC+=RxData[index+i];
					}
					if(SUMCRC!=RxData[index+10]){
						index+=11;						
						continue;
					}
					else{
						raw_Roll=(short)((short)RxData[index+3]<<8|RxData[index+2]);
						raw_Pitch=(short)((short)RxData[index+5]<<8|RxData[index+4]);
						raw_Yaw=(short)((short)RxData[index+7]<<8|RxData[index+6]);
					}					
					break;
				case 0x59:
					for(int i=0;i<10;i++){
					SUMCRC+=RxData[index+i];
					}
					if(SUMCRC!=RxData[index+10]){
						index+=11;						
						continue;
					}
					else{
						raw_w=(short)((short)RxData[index+3]<<8|RxData[index+2]);
						raw_x=(short)((short)RxData[index+5]<<8|RxData[index+4]);
						raw_y=(short)((short)RxData[index+7]<<8|RxData[index+6]);
						raw_z=(short)((short)RxData[index+9]<<8|RxData[index+8]);
					}					
					break;					
				default:
					index+=11;
					continue;
			}
			index+=11;
		}
//		printf("ax:%hd\tay:%hd\taz:%hd\r\n",ax,ay,az);
//		printf("wx:%hd\twy:%hd\twz:%hd\r\n",wx,wy,wz);
//		printf("ROLL:%hd\tPITCH:%hd\tYAW:%hd\r\n",Roll,Pitch,Yaw);
		real_ax=(1.0*raw_ax*16*g)/32768;
		real_ay=(1.0*raw_ay*16*g)/32768;
		real_az=(1.0*raw_az*16*g)/32768;
		real_wx=(1.0*raw_wx*2000)/32768;
		real_wy=(1.0*raw_wy*2000)/32768;
		real_wz=(1.0*raw_wz*2000)/32768;
		real_Roll=(1.0*raw_Roll*180)/32768;
		real_Pitch=(1.0*raw_Pitch*180)/32768;
		real_Yaw=(1.0*raw_Yaw*180)/32768;
		real_w=1.0*raw_w/32768;
		real_x=1.0*raw_x/32768;
		real_y=1.0*raw_y/32768;
		real_z=1.0*raw_z/32768;
		
		convert_body_to_world_acceleration(real_ax,real_ay,real_az,real_w,real_x,real_y,real_z,&world_ax,&world_ay,&world_az);
		if(zero_offset_flat<100)
		{
			a_bias[0] += world_ax;
			a_bias[1] += world_ay;
			a_bias[2] += world_az - 9.8f; // 假设Z轴向下，减去重力
			zero_offset_flat++;
			if(zero_offset_flat==100)
			{
				a_bias[0] /= 100; a_bias[1] /= 100; a_bias[2] /= 100;
				zero_offset_flat++;
			}
		}

		if(zero_offset_flat>=100)
		{
		float a_corrected[3];
//		a_corrected[0] = world_ax - a_bias[0];
//		a_corrected[1] = world_ay - a_bias[1];
//		a_corrected[2] = world_az - a_bias[2] - 9.8f; // 再次确保去除重力

		a_corrected[0] = world_ax;
		a_corrected[1] = world_ay;
		a_corrected[2] = world_az - 9.8f; // 再次确保去除重力
		// 积分得到速度
//		filtered_ax = alpha * (filtered_ax + a_corrected[0] * DT) - alpha * a_corrected[0];
//        filtered_ay = alpha * (filtered_ay + a_corrected[1] * DT) - alpha * a_corrected[1];
//        filtered_az = alpha * (filtered_az + a_corrected[2] * DT) - alpha * a_corrected[2];

//		float acc_magnitude = sqrt(filtered_ax * filtered_ax + filtered_ay * filtered_ay + filtered_az * filtered_az);
		float acc_magnitude = sqrt(a_corrected[0]*a_corrected[0] + a_corrected[1]*a_corrected[1] + a_corrected[2]*a_corrected[2]);
        float gyro_magnitude = sqrt(real_wx * real_wx + real_wy * real_wy + real_wz * real_wz);

		static float prev_vz = 0;
		int direction_changed = 0;
    
    // 预测当前速度变化
    float predicted_vz = v_prev[2] + a_corrected[2] * DT;
    if(predicted_vz * prev_vz < 0 && fabs(predicted_vz) > 0.1f) {
        direction_changed = 1;
    }
    prev_vz = predicted_vz;

        // 严格的静止检测阈值
        if((acc_magnitude < 1.0f && gyro_magnitude < 20.0f)||direction_changed) {
            stop_cnt++;
        } else {
            stop_cnt = 0;
        }

        float v_current[3];
        if(stop_cnt >= 3) { // 连续5帧静止才重置速度
            v_current[0] = 0;
            v_current[1] = 0;
            v_current[2] = 0;
			
//			filtered_ax = 0;
//			filtered_ay = 0;
//			filtered_az = 0;
            // 可选：在静止时轻微调整零偏
//            a_bias[0] = 0.999f * a_bias[0] + 0.001f * world_ax;
//            a_bias[1] = 0.999f * a_bias[1] + 0.001f * world_ay;
//            a_bias[2] = 0.999f * a_bias[2] + 0.001f * (world_az + 9.8f); // 对应校准公式
        } else {
            // 运动状态，使用滤波后的加速度进行积分
		v_current[0] = v_prev[0] + a_corrected[0] * DT;
		v_current[1] = v_prev[1] + a_corrected[1] * DT;
		v_current[2] = v_prev[2] + a_corrected[2] * DT;
//		v_current[0] = v_prev[0] + filtered_ax * DT;
//        v_current[1] = v_prev[1] + filtered_ay * DT;
//        v_current[2] = v_prev[2] + filtered_az * DT;
        }
		// 积分得到位移
		float s_current[3];			
		s_current[0] = s_prev[0] + v_current[0] * DT;
		s_current[1] = s_prev[1] + v_current[1] * DT;
		s_current[2] = s_prev[2] + v_current[2] * DT;

		// 更新状态
		v_prev[0] = v_current[0];
		v_prev[1] = v_current[1];
		v_prev[2] = v_current[2];
		s_prev[0] = s_current[0];
		s_prev[1] = s_current[1];
		s_prev[2] = s_current[2];
		



	//		printf("ax:%f\tay:%f\taz:%f\twx:%f\twy:%f\twz:%f\tROLL:%f\tPITCH:%f\tYAW:%f\r\n",real_ax
	//		,real_ay,real_az,real_wx,real_wy,real_wz,real_Roll,real_Pitch,real_Yaw);
			sendSensorData(a_corrected[0],a_corrected[1],a_corrected[2],v_prev[0],v_prev[1],v_prev[2],s_prev[0],s_prev[1],s_prev[2],
			real_wx,real_wy,real_wz,direction_changed);
	//		sendSensorData(real_ax,real_ay,real_az,real_wx,real_wy,real_wz,real_Roll,real_Pitch,real_Yaw,
	//		real_w,real_x,real_y,real_z);
	//		float test1=1;
	//		sendSensorData(test1,test1,test1,test1,test1,test1,test1,test1,test1,test1,test1,test1,test1);

		}
}

