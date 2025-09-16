#include "matrix.h"
#include <string.h> // 用于memcpy

/**
 * @brief 初始化矩阵（不填充特定值）
 * @param mat 指向矩阵的指针
 * @param rows 矩阵行数
 * @param cols 矩阵列数
 */
void matrix_init(Matrix* mat, uint8_t rows, uint8_t cols) {
    // 检查并限制矩阵尺寸不超过最大值
    if (rows > MAX_MATRIX_ROWS) rows = MAX_MATRIX_ROWS;
    if (cols > MAX_MATRIX_COLS) cols = MAX_MATRIX_COLS;
    
    mat->rows = rows;
    mat->cols = cols;
    
    // 将所有元素初始化为0
    for (uint8_t i = 0; i < rows; i++) {
        for (uint8_t j = 0; j < cols; j++) {
            mat->data[i][j] = 0.0f;
        }
    }
}

/**
 * @brief 初始化零矩阵（所有元素设为0）
 * @param mat 指向矩阵的指针
 * @param rows 矩阵行数
 * @param cols 矩阵列数
 */
void matrix_init_zeros(Matrix* mat, uint8_t rows, uint8_t cols) {
    matrix_init(mat, rows, cols); // 默认就是零矩阵
}

/**
 * @brief 初始化单位矩阵（对角线为1，其余为0）
 * @param mat 指向矩阵的指针
 * @param size 矩阵尺寸（方阵）
 */
void matrix_init_identity(Matrix* mat, uint8_t size) {
    // 检查并限制矩阵尺寸不超过最大值
    if (size > MAX_MATRIX_ROWS) size = MAX_MATRIX_ROWS;
    if (size > MAX_MATRIX_COLS) size = MAX_MATRIX_COLS;
    
    mat->rows = size;
    mat->cols = size;
    
    // 将对角线元素设为1，其余设为0
    for (uint8_t i = 0; i < size; i++) {
        for (uint8_t j = 0; j < size; j++) {
            mat->data[i][j] = (i == j) ? 1.0f : 0.0f;
        }
    }
}

/**
 * @brief 复制矩阵
 * @param src 源矩阵指针
 * @param dest 目标矩阵指针
 */
void matrix_copy(const Matrix* src, Matrix* dest) {
    dest->rows = src->rows;
    dest->cols = src->cols;
    
    // 复制所有元素
    for (uint8_t i = 0; i < src->rows; i++) {
        for (uint8_t j = 0; j < src->cols; j++) {
            dest->data[i][j] = src->data[i][j];
        }
    }
}

/**
 * @brief 打印矩阵内容（需要通过串口实现）
 * @param mat 要打印的矩阵指针
 */
void matrix_print(const Matrix* mat) {
    // 需要根据您的串口实现来编写
    // 例如: printf("Matrix %dx%d:\n", mat->rows, mat->cols);
    // for (uint8_t i = 0; i < mat->rows; i++) {
    //     for (uint8_t j = 0; j < mat->cols; j++) {
    //         printf("%.4f\t", mat->data[i][j]);
    //     }
    //     printf("\n");
    // }
}

/**
 * @brief 矩阵加法运算
 * @param a 第一个矩阵指针
 * @param b 第二个矩阵指针
 * @param res 结果矩阵指针
 * @return 0表示成功，1表示维度不匹配
 */
uint8_t matrix_add(const Matrix* a, const Matrix* b, Matrix* res) {
    // 检查矩阵维度是否匹配
    if (a->rows != b->rows || a->cols != b->cols) {
        return 1; // 维度不匹配
    }
    
    res->rows = a->rows;
    res->cols = a->cols;
    
    // 执行逐元素加法
    for (uint8_t i = 0; i < a->rows; i++) {
        for (uint8_t j = 0; j < a->cols; j++) {
            res->data[i][j] = a->data[i][j] + b->data[i][j];
        }
    }
    
    return 0; // 成功
}

/**
 * @brief 矩阵减法运算
 * @param a 第一个矩阵指针
 * @param b 第二个矩阵指针
 * @param res 结果矩阵指针
 * @return 0表示成功，1表示维度不匹配
 */
uint8_t matrix_subtract(const Matrix* a, const Matrix* b, Matrix* res) {
    // 检查矩阵维度是否匹配
    if (a->rows != b->rows || a->cols != b->cols) {
        return 1; // 维度不匹配
    }
    
    res->rows = a->rows;
    res->cols = a->cols;
    
    // 执行逐元素减法
    for (uint8_t i = 0; i < a->rows; i++) {
        for (uint8_t j = 0; j < a->cols; j++) {
            res->data[i][j] = a->data[i][j] - b->data[i][j];
        }
    }
    
    return 0; // 成功
}

/**
 * @brief 矩阵乘法运算
 * @param a 第一个矩阵指针
 * @param b 第二个矩阵指针
 * @param res 结果矩阵指针
 * @return 0表示成功，1表示维度不匹配
 */
uint8_t matrix_multiply(const Matrix* a, const Matrix* b, Matrix* res) {
    // 检查矩阵维度是否匹配（a的列数必须等于b的行数）
    if (a->cols != b->rows) {
        return 1; // 维度不匹配
    }
    
    res->rows = a->rows;
    res->cols = b->cols;
    
    // 执行矩阵乘法
    for (uint8_t i = 0; i < a->rows; i++) {
        for (uint8_t j = 0; j < b->cols; j++) {
            float sum = 0.0f;
            for (uint8_t k = 0; k < a->cols; k++) {
                sum += a->data[i][k] * b->data[k][j];
            }
            res->data[i][j] = sum;
        }
    }
    
    return 0; // 成功
}

/**
 * @brief 矩阵标量乘法
 * @param mat 矩阵指针
 * @param scalar 标量值
 * @param res 结果矩阵指针
 * @return 0表示成功
 */
uint8_t matrix_multiply_scalar(const Matrix* mat, float scalar, Matrix* res) {
    res->rows = mat->rows;
    res->cols = mat->cols;
    
    // 执行标量乘法
    for (uint8_t i = 0; i < mat->rows; i++) {
        for (uint8_t j = 0; j < mat->cols; j++) {
            res->data[i][j] = mat->data[i][j] * scalar;
        }
    }
    
    return 0; // 成功
}

/**
 * @brief 矩阵转置
 * @param mat 输入矩阵指针
 * @param res 结果矩阵指针
 * @return 0表示成功
 */
uint8_t matrix_transpose(const Matrix* mat, Matrix* res) {
    res->rows = mat->cols;
    res->cols = mat->rows;
    
    // 执行转置操作
    for (uint8_t i = 0; i < mat->rows; i++) {
        for (uint8_t j = 0; j < mat->cols; j++) {
            res->data[j][i] = mat->data[i][j];
        }
    }
    
    return 0; // 成功
}

/**
 * @brief 3x3矩阵求逆
 * @param mat 输入矩阵指针（必须是3x3）
 * @param res 结果矩阵指针
 * @return 0表示成功，1表示不是3x3矩阵，2表示矩阵不可逆
 */
uint8_t matrix_3x3_inverse(const Matrix* mat, Matrix* res) {
    // 检查矩阵是否为3x3
    if (mat->rows != 3 || mat->cols != 3) {
        return 1; // 不是3x3矩阵
    }
    
    // 计算行列式
    float det = matrix_3x3_determinant(mat);
    
    // 检查行列式是否接近零（矩阵是否可逆）
    if (fabsf(det) < 1e-6f) {
        return 2; // 矩阵不可逆
    }
    
    float inv_det = 1.0f / det;
    
    // 计算伴随矩阵并转置（即逆矩阵）
    res->data[0][0] = (mat->data[1][1] * mat->data[2][2] - mat->data[1][2] * mat->data[2][1]) * inv_det;
    res->data[0][1] = (mat->data[0][2] * mat->data[2][1] - mat->data[0][1] * mat->data[2][2]) * inv_det;
    res->data[0][2] = (mat->data[0][1] * mat->data[1][2] - mat->data[0][2] * mat->data[1][1]) * inv_det;
    
    res->data[1][0] = (mat->data[1][2] * mat->data[2][0] - mat->data[1][0] * mat->data[2][2]) * inv_det;
    res->data[1][1] = (mat->data[0][0] * mat->data[2][2] - mat->data[0][2] * mat->data[2][0]) * inv_det;
    res->data[1][2] = (mat->data[0][2] * mat->data[1][0] - mat->data[0][0] * mat->data[1][2]) * inv_det;
    
    res->data[2][0] = (mat->data[1][0] * mat->data[2][1] - mat->data[1][1] * mat->data[2][0]) * inv_det;
    res->data[2][1] = (mat->data[0][1] * mat->data[2][0] - mat->data[0][0] * mat->data[2][1]) * inv_det;
    res->data[2][2] = (mat->data[0][0] * mat->data[1][1] - mat->data[0][1] * mat->data[1][0]) * inv_det;
    
    res->rows = 3;
    res->cols = 3;
    
    return 0; // 成功
}

/**
 * @brief 计算3x3矩阵的行列式
 * @param mat 输入矩阵指针（必须是3x3）
 * @return 行列式值，如果不是3x3矩阵则返回0
 */
float matrix_3x3_determinant(const Matrix* mat) {
    // 检查矩阵是否为3x3
    if (mat->rows != 3 || mat->cols != 3) {
        return 0.0f; // 不是3x3矩阵
    }
    
    // 计算3x3矩阵的行列式
    return mat->data[0][0] * (mat->data[1][1] * mat->data[2][2] - mat->data[1][2] * mat->data[2][1]) -
           mat->data[0][1] * (mat->data[1][0] * mat->data[2][2] - mat->data[1][2] * mat->data[2][0]) +
           mat->data[0][2] * (mat->data[1][0] * mat->data[2][1] - mat->data[1][1] * mat->data[2][0]);
}

/**
 * @brief 比较两个矩阵是否相等（在容差范围内）
 * @param a 第一个矩阵指针
 * @param b 第二个矩阵指针
 * @param tolerance 容差值
 * @return 1表示相等，0表示不相等
 */
uint8_t matrix_equal(const Matrix* a, const Matrix* b, float tolerance) {
    // 检查矩阵维度是否相同
    if (a->rows != b->rows || a->cols != b->cols) {
        return 0; // 维度不同
    }
    
    // 比较所有元素是否在容差范围内相等
    for (uint8_t i = 0; i < a->rows; i++) {
        for (uint8_t j = 0; j < a->cols; j++) {
            if (fabsf(a->data[i][j] - b->data[i][j]) > tolerance) {
                return 0; // 元素不相等
            }
        }
    }
    
    return 1; // 矩阵相等
}

/**
 * @brief 设置矩阵元素值
 * @param mat 矩阵指针
 * @param row 行索引
 * @param col 列索引
 * @param value 要设置的值
 */
void matrix_set(Matrix* mat, uint8_t row, uint8_t col, float value) {
    // 检查索引是否有效
    if (row < mat->rows && col < mat->cols) {
        mat->data[row][col] = value;
    }
}

/**
 * @brief 获取矩阵元素值
 * @param mat 矩阵指针
 * @param row 行索引
 * @param col 列索引
 * @return 元素值，如果索引越界则返回0
 */
float matrix_get(const Matrix* mat, uint8_t row, uint8_t col) {
    // 检查索引是否有效
    if (row < mat->rows && col < mat->cols) {
        return mat->data[row][col];
    }
    return 0.0f;
}

/**
 * @brief 将四元数转换为旋转矩阵
 * @param qw 四元数w分量
 * @param qx 四元数x分量
 * @param qy 四元数y分量
 * @param qz 四元数z分量
 * @param rot_mat 旋转矩阵指针（3x3）
 */
void quaternion_to_rotation_matrix(float qw, float qx, float qy, float qz, Matrix* rot_mat) {
    // 确保四元数是单位四元数
    float norm = sqrtf(qw*qw + qx*qx + qy*qy + qz*qz);
    if (fabsf(norm - 1.0f) > 1e-6f) {
        qw /= norm;
        qx /= norm;
        qy /= norm;
        qz /= norm;
    }
    
    // 计算旋转矩阵元素
    rot_mat->data[0][0] = 1.0f - 2.0f*(qy*qy + qz*qz);
    rot_mat->data[0][1] = 2.0f*(qx*qy - qw*qz);
    rot_mat->data[0][2] = 2.0f*(qx*qz + qw*qy);
    
    rot_mat->data[1][0] = 2.0f*(qx*qy + qw*qz);
    rot_mat->data[1][1] = 1.0f - 2.0f*(qx*qx + qz*qz);
    rot_mat->data[1][2] = 2.0f*(qy*qz - qw*qx);
    
    rot_mat->data[2][0] = 2.0f*(qx*qz - qw*qy);
    rot_mat->data[2][1] = 2.0f*(qy*qz + qw*qx);
    rot_mat->data[2][2] = 1.0f - 2.0f*(qx*qx + qy*qy);
    
    rot_mat->rows = 3;
    rot_mat->cols = 3;
}

/**
 * @brief 使用四元数旋转向量
 * @param vec 输入向量（3元素数组）
 * @param qw 四元数w分量
 * @param qx 四元数x分量
 * @param qy 四元数y分量
 * @param qz 四元数z分量
 * @param result 结果向量（3元素数组）
 */
void rotate_vector_by_quaternion(const float vec[3], float qw, float qx, float qy, float qz, float result[3]) {
    // 确保四元数是单位四元数
    float norm = sqrtf(qw*qw + qx*qx + qy*qy + qz*qz);
    if (fabsf(norm - 1.0f) > 1e-6f) {
        qw /= norm;
        qx /= norm;
        qy /= norm;
        qz /= norm;
    }
    
    // 直接计算旋转后的向量（优化算法）
    float ux = qx, uy = qy, uz = qz;
    float s = qw;
    
    // 计算 cross(u, v)
    float cross_uv_x = uy * vec[2] - uz * vec[1];
    float cross_uv_y = uz * vec[0] - ux * vec[2];
    float cross_uv_z = ux * vec[1] - uy * vec[0];
    
    // 计算 cross(u, cross(u, v)) + s * cross(u, v)
    float cross_u_cross_uv_x = uy * cross_uv_z - uz * cross_uv_y;
    float cross_u_cross_uv_y = uz * cross_uv_x - ux * cross_uv_z;
    float cross_u_cross_uv_z = ux * cross_uv_y - uy * cross_uv_x;
    
    result[0] = vec[0] + 2.0f * (cross_u_cross_uv_x + s * cross_uv_x);
    result[1] = vec[1] + 2.0f * (cross_u_cross_uv_y + s * cross_uv_y);
    result[2] = vec[2] + 2.0f * (cross_u_cross_uv_z + s * cross_uv_z);
}
