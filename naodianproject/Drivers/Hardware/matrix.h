#ifndef MATRIX_H
#define MATRIX_H
#include "stm32f4xx.h"                  // Device header
#include <math.h>

// 定义最大矩阵尺寸（根据您的需求调整）
#define MAX_MATRIX_ROWS 4
#define MAX_MATRIX_COLS 4

// 矩阵结构体
typedef struct {
    uint8_t rows;
    uint8_t cols;
    float data[MAX_MATRIX_ROWS][MAX_MATRIX_COLS];
} Matrix;

// 初始化函数
void matrix_init(Matrix* mat, uint8_t rows, uint8_t cols);
void matrix_init_zeros(Matrix* mat, uint8_t rows, uint8_t cols);
void matrix_init_identity(Matrix* mat, uint8_t size);

// 基本操作
void matrix_copy(const Matrix* src, Matrix* dest);
void matrix_print(const Matrix* mat); // 需要通过串口实现

// 矩阵运算 (结果存储在res中，返回0表示成功，非0表示错误)
uint8_t matrix_add(const Matrix* a, const Matrix* b, Matrix* res);
uint8_t matrix_subtract(const Matrix* a, const Matrix* b, Matrix* res);
uint8_t matrix_multiply(const Matrix* a, const Matrix* b, Matrix* res);
uint8_t matrix_multiply_scalar(const Matrix* mat, float scalar, Matrix* res);
uint8_t matrix_transpose(const Matrix* mat, Matrix* res);

// 3x3矩阵专用函数（更高效）
uint8_t matrix_3x3_inverse(const Matrix* mat, Matrix* res);
float matrix_3x3_determinant(const Matrix* mat);

// 实用函数
uint8_t matrix_equal(const Matrix* a, const Matrix* b, float tolerance);
void matrix_set(Matrix* mat, uint8_t row, uint8_t col, float value);
float matrix_get(const Matrix* mat, uint8_t row, uint8_t col);

// 四元数相关函数
void quaternion_to_rotation_matrix(float qw, float qx, float qy, float qz, Matrix* rot_mat);
void rotate_vector_by_quaternion(const float vec[3], float qw, float qx, float qy, float qz, float result[3]);

#endif
