#include "datastructure.h"

USHORT RGB16Bit565(int r, int g, int b)
{
    // this function simply builds a 5.6.5 format 16 bit pixel
    // assumes input is RGB 0-255 each channel
    r >>= 3;
    g >>= 2;
    b >>= 3;
    return (_RGB16BIT565((r), (g), (b)));
}

float *cos_look;
float *sin_look;

void Build_Sin_Cos_Tables(void)
{
    cos_look = new float[361];
    sin_look = new float[361];
    for (int ang = 0; ang <= 360; ang++)
    {
        float theta = (float)ang * PI / (float)180;
        cos_look[ang] = cos(theta);
        sin_look[ang] = sin(theta);
    }
}

float Fast_Sin(float theta)
{
    // this function uses the sin_look[] lookup table, but
    // has logic to handle negative angles as well as fractional
    // angles via interpolation, use this for a more robust
    // sin computation that the blind lookup, but with with
    // a slight hit in speed

    // convert angle to 0-359
    theta = fmodf(theta, 360);

    // make angle positive
    if (theta < 0)
        theta += 360.0;

    // compute floor of theta and fractional part to interpolate
    int theta_int = (int)theta;
    float theta_frac = theta - theta_int;

    // now compute the value of sin(angle) using the lookup tables
    // and interpolating the fractional part, note that if theta_int
    // is equal to 359 then theta_int+1=360, but this is fine since the
    // table was made with the entries 0-360 inclusive
    return (sin_look[theta_int] + theta_frac * (sin_look[theta_int + 1] - sin_look[theta_int]));

} // end Fast_Sin

///////////////////////////////////////////////////////////////

float Fast_Cos(float theta)
{
    // this function uses the cos_look[] lookup table, but
    // has logic to handle negative angles as well as fractional
    // angles via interpolation, use this for a more robust
    // cos computation that the blind lookup, but with with
    // a slight hit in speed

    // convert angle to 0-359
    theta = fmodf(theta, 360);

    // make angle positive
    if (theta < 0)
        theta += 360.0;

    // compute floor of theta and fractional part to interpolate
    int theta_int = (int)theta;
    float theta_frac = theta - theta_int;

    // now compute the value of sin(angle) using the lookup tables
    // and interpolating the fractional part, note that if theta_int
    // is equal to 359 then theta_int+1=360, but this is fine since the
    // table was made with the entries 0-360 inclusive
    return (cos_look[theta_int] + theta_frac * (cos_look[theta_int + 1] - cos_look[theta_int]));

} // end Fast_Cos

void PLANE3D_Init(PLANE3D_PTR plane, POINT3D_PTR p0,
                  VECTOR3D_PTR normal, int normalize = 0)
{
    // this function initializes a 3d plane

    // copy the point
    POINT3D_COPY(&plane->p0, p0);

    // if normalize is 1 then the normal is made into a unit vector
    if (!normalize)
        VECTOR3D_COPY(&plane->n, normal);
    else
    {
        // make normal into unit vector
        VECTOR3D_Normalize(normal, &plane->n);
    } // end else
}

void VECTOR3D_Normalize(VECTOR3D_PTR va)
{
    // normalizes the sent vector in placew

    // compute length
    float length = sqrtf(va->x * va->x + va->y * va->y + va->z * va->z);

    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;

    float length_inv = 1 / length;

    // compute normalized version of vector
    va->x *= length_inv;
    va->y *= length_inv;
    va->z *= length_inv;

} // end VECTOR3D_Normalize

//向量正则化
void VECTOR3D_Normalize(VECTOR3D_PTR va, VECTOR3D_PTR vn)
{
    // normalizes the sent vector and returns the result in vn

    VECTOR3D_ZERO(vn);

    // compute length
    float length = VECTOR3D_Length(va);

    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;

    float length_inv = 1.0 / length;

    // compute normalized version of vector
    vn->x = va->x * length_inv;
    vn->y = va->y * length_inv;
    vn->z = va->z * length_inv;

} // end VECTOR3D_Normalize

//向量长度
float VECTOR3D_Length(VECTOR3D_PTR va)
{
    // computes the magnitude of a vector, slow

    return ((float)sqrtf(va->x * va->x + va->y * va->y + va->z * va->z));

} // end VECTOR3D_Length

//构造旋转矩阵，前3个参数是三个轴的旋转角度，最后一个参数是旋转矩阵的矩阵指针
void Build_XYZ_Rotation_MATRIX4X4(float theta_x, float theta_y, float theta_z, MATRIX4X4_PTR mrot)
{
    MATRIX4X4 mx, my, mz, mtmp;         // working matrices
    float sin_theta = 0, cos_theta = 0; // used to initialize matrices
    int rot_seq = 0;                    // 1 for x, 2 for y, 4 for z

    // step 0: fill in with identity matrix
    MAT_IDENTITY_4X4(mrot);

    // step 1: based on zero and non-zero rotation angles, determine
    // rotation sequence
    if (fabs(theta_x) > EPSILON_E5) // x
        rot_seq = rot_seq | 1;      //0001

    if (fabs(theta_y) > EPSILON_E5) // y
        rot_seq = rot_seq | 2;      //0010

    if (fabs(theta_z) > EPSILON_E5) // z
        rot_seq = rot_seq | 4;      //0100

    switch (rot_seq)
    {
    case 0:
    {
        return;
    }
    break;

    case 1: // x rotation
    {
        // compute the sine and cosine of the angle
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        MAT_COPY_4X4(&mx, mrot);
        return;
    }
    break;

    case 2: // y rotation
    {
        // compute the sine and cosine of the angle
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my,
                     cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        MAT_COPY_4X4(&my, mrot);
        return;
    }
    break;

    case 3: // xy rotation
    {
        // compute the sine and cosine of the angle for x
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle for y
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // concatenate matrices
        Mat_Mul_4X4(&mx, &my, mrot);
        return;
    }
    break;

    case 4: // z rotation
    {
        // compute the sine and cosine of the angle
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // that's it, copy to output matrix
        MAT_COPY_4X4(&mz, mrot);
        return;
    }
    break;

    case 5: // xz rotation
    {
        // compute the sine and cosine of the angle x
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle z
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // concatenate matrices
        Mat_Mul_4X4(&mx, &mz, mrot);
        return;
    }
    break;

    case 6: // yz rotation
    {
        // compute the sine and cosine of the angle y
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle z
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // concatenate matrices
        Mat_Mul_4X4(&my, &mz, mrot);
        return;
    }
    break;

    case 7: // xyz rotation
    {
        // compute the sine and cosine of the angle x
        cos_theta = Fast_Cos(theta_x);
        sin_theta = Fast_Sin(theta_x);

        // set the matrix up
        Mat_Init_4X4(&mx, 1, 0, 0, 0,
                     0, cos_theta, sin_theta, 0,
                     0, -sin_theta, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle y
        cos_theta = Fast_Cos(theta_y);
        sin_theta = Fast_Sin(theta_y);

        // set the matrix up
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // compute the sine and cosine of the angle z
        cos_theta = Fast_Cos(theta_z);
        sin_theta = Fast_Sin(theta_z);

        // set the matrix up
        Mat_Init_4X4(&mz, cos_theta, sin_theta, 0, 0,
                     -sin_theta, cos_theta, 0, 0,
                     0, 0, 1, 0,
                     0, 0, 0, 1);

        // concatenate matrices, watch order!
        Mat_Mul_4X4(&mx, &my, &mtmp);
        Mat_Mul_4X4(&mtmp, &mz, mrot);
    }
    break;

    default:
        break;

    } // end switch

} // end Build_XYZ_Rotation_MATRIX4X4

void MAT_IDENTITY_4X4(MATRIX4X4_PTR m)
{
    memcpy((void *)(m), (void *)&IMAT_4X4, sizeof(MATRIX4X4));
}

//向量va 乘 矩阵mv，得到新的向量vprod
void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR va, MATRIX4X4_PTR mb, VECTOR4D_PTR vprod)
{
    for (int col = 0; col < 4; col++)
    {
        float sum = 0;
        for (int row = 0; row < 4; row++)
        {
            sum += (va->M[row] * mb->M[row][col]);
        }
        vprod->M[col] = sum;
    }
}

//向量加法
//vsum = va + vb
void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum)
{
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
    vsum->w = 1;
}

//向量加法
VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    VECTOR4D vsum;

    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;
    vsum.w = 1;
    return (vsum);
}

//矩阵相乘
void Mat_Mul_4X4(MATRIX4X4_PTR ma, MATRIX4X4_PTR mb, MATRIX4X4_PTR mprod)
{
    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 4; col++)
        {
            float sum = 0;
            for (int index = 0; index < 4; index++)
            {
                sum += (ma->M[row][index] * mb->M[index][col]);
            }
            mprod->M[row][col] = sum;
        }
    }
}

//矩阵初始化
void Mat_Init_4X4(MATRIX4X4_PTR ma,
                  float m00, float m01, float m02, float m03,
                  float m10, float m11, float m12, float m13,
                  float m20, float m21, float m22, float m23,
                  float m30, float m31, float m32, float m33)

{
    // this function fills a 4x4 matrix with the sent data in
    // row major form
    ma->M00 = m00;
    ma->M01 = m01;
    ma->M02 = m02;
    ma->M03 = m03;
    ma->M10 = m10;
    ma->M11 = m11;
    ma->M12 = m12;
    ma->M13 = m13;
    ma->M20 = m20;
    ma->M21 = m21;
    ma->M22 = m22;
    ma->M23 = m23;
    ma->M30 = m30;
    ma->M31 = m31;
    ma->M32 = m32;
    ma->M33 = m33;

} // end Mat_Init_4X4

void Init_CAM4DV1(CAM4DV1_PTR cam,
                  int cam_attr,
                  POINT4D_PTR cam_pos,    // 相机初始位置
                  VECTOR4D_PTR cam_dir,   // 相机初始角度
                  POINT4D_PTR cam_target, // UVN相机的目标位置
                  float near_clip_z,      // 近剪裁面和远剪裁面
                  float far_clip_z,
                  float fov,            // 视野，单位为度
                  float viewport_width, // 屏幕宽高
                  float viewport_height)
{

    cam->attr = cam_attr;

    VECTOR4D_COPY(&cam->pos, cam_pos);
    VECTOR4D_COPY(&cam->dir, cam_dir);

    // for UVN camera
    VECTOR4D_INITXYZ(&cam->u, 1, 0, 0); // set to +x
    VECTOR4D_INITXYZ(&cam->v, 0, 1, 0); // set to +y
    VECTOR4D_INITXYZ(&cam->n, 0, 0, 1); // set to +z

    if (cam_target != nullptr)
        VECTOR4D_COPY(&cam->target, cam_target); // UVN target
    else
        VECTOR4D_ZERO(&cam->target);

    cam->near_clip_z = near_clip_z; // 近剪裁面
    cam->far_clip_z = far_clip_z;   // 远剪裁面

    cam->viewport_width = viewport_width;
    cam->viewport_height = viewport_height;

    cam->viewport_center_x = (viewport_width - 1) / 2;
    cam->viewport_center_y = (viewport_height - 1) / 2;
    cam->aspect_ratio = (float)viewport_width / (float)viewport_height;

    // 将所有矩阵设置为单位矩阵
    MAT_IDENTITY_4X4(&cam->mcam);
    MAT_IDENTITY_4X4(&cam->mper);
    MAT_IDENTITY_4X4(&cam->mscr);

    cam->fov = fov;

    // 设置视平面，高为2 x (2/对比度)
    cam->viewplane_width = 2.0;
    cam->viewplane_height = 2.0 / cam->aspect_ratio;

    // now we know fov and we know the viewplane dimensions plug into formula and
    // solve for view distance parameters
    float tan_fov_div2 = tan(DEG_TO_RAD(fov / 2));

    cam->view_dist = (0.5) * (cam->viewplane_width) * tan_fov_div2;

    // test for 90 fov first since it's easy :)
    if (fov == 90.0)
    {
        // set up the clipping planes -- easy for 90 degrees!
        POINT3D pt_origin; // point on the plane
        VECTOR3D_INITXYZ(&pt_origin, 0, 0, 0);

        VECTOR3D vn; // normal to plane

        // right clipping plane
        VECTOR3D_INITXYZ(&vn, 1, 0, -1); // x=z plane
        PLANE3D_Init(&cam->rt_clip_plane, &pt_origin, &vn, 1);

        // left clipping plane
        VECTOR3D_INITXYZ(&vn, -1, 0, -1); // -x=z plane
        PLANE3D_Init(&cam->lt_clip_plane, &pt_origin, &vn, 1);

        // top clipping plane
        VECTOR3D_INITXYZ(&vn, 0, 1, -1); // y=z plane
        PLANE3D_Init(&cam->tp_clip_plane, &pt_origin, &vn, 1);

        // bottom clipping plane
        VECTOR3D_INITXYZ(&vn, 0, -1, -1); // -y=z plane
        PLANE3D_Init(&cam->bt_clip_plane, &pt_origin, &vn, 1);
    } // end if d=1
    else
    {
        // now compute clipping planes yuck!
        POINT3D pt_origin; // point on the plane
        VECTOR3D_INITXYZ(&pt_origin, 0, 0, 0);

        VECTOR3D vn; // normal to plane

        // since we don't have a 90 fov, computing the normals
        // are a bit tricky, there are a number of geometric constructions
        // that solve the problem, but I'm going to solve for the
        // vectors that represent the 2D projections of the frustrum planes
        // on the x-z and y-z planes and then find perpendiculars to them

        // right clipping plane, check the math on graph paper
        VECTOR3D_INITXYZ(&vn, cam->view_dist, 0, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->rt_clip_plane, &pt_origin, &vn, 1);

        // left clipping plane, we can simply reflect the right normal about
        // the z axis since the planes are symetric about the z axis
        // thus invert x only
        VECTOR3D_INITXYZ(&vn, -cam->view_dist, 0, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->lt_clip_plane, &pt_origin, &vn, 1);

        // top clipping plane, same construction
        VECTOR3D_INITXYZ(&vn, 0, cam->view_dist, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->tp_clip_plane, &pt_origin, &vn, 1);

        // bottom clipping plane, same inversion
        VECTOR3D_INITXYZ(&vn, 0, -cam->view_dist, -cam->viewplane_width / 2.0);
        PLANE3D_Init(&cam->bt_clip_plane, &pt_origin, &vn, 1);
    } // end else

} // end Init_CAM4DV1

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                    POLYF4DV1_PTR poly)
{
    if (rend_list->num_polys >= RENDERLIST4DV1_MAX_POLYS)
        return (0);

    rend_list->poly_ptrs[rend_list->num_polys] = &rend_list->poly_data[rend_list->num_polys];

    memcpy((void *)&rend_list->poly_data[rend_list->num_polys], (void *)poly, sizeof(POLYF4DV1));

    if (rend_list->num_polys == 0)
    {
        rend_list->poly_data[0].next = NULL;
        rend_list->poly_data[0].prev = NULL;
    }
    else
    {
        rend_list->poly_data[rend_list->num_polys].next = NULL;
        rend_list->poly_data[rend_list->num_polys].prev =
            &rend_list->poly_data[rend_list->num_polys - 1];

        rend_list->poly_data[rend_list->num_polys - 1].next =
            &rend_list->poly_data[rend_list->num_polys];
    }
    rend_list->num_polys++;
    return (1);
}

//构造欧拉相机矩阵
void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq)
{
    MATRIX4X4 mt_inv, // 平移矩阵
        mx_inv,       //  绕x轴的旋转矩阵
        my_inv,       //  绕y轴的旋转矩阵
        mz_inv,       //  绕z轴的旋转矩阵
        mrot,         //  总的旋转矩阵
        mtmp;         //    临时矩阵

    // 第一步：平移矩阵
    Mat_Init_4X4(&mt_inv,
                 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 -cam->pos.x, -cam->pos.y, -cam->pos.z, 1);

    // 第二步：三个旋转矩阵
    float theta_x = cam->dir.x;
    float theta_y = cam->dir.y;
    float theta_z = cam->dir.z;

    //  绕x轴的旋转矩阵
    float cos_theta = Fast_Cos(theta_x);        // cos(-x) = cos(x)
    float minus_sin_theta = -Fast_Sin(theta_x); // sin(-x) = -sin(x)
    Mat_Init_4X4(&mx_inv,
                 1, 0, 0, 0,
                 0, cos_theta, minus_sin_theta, 0,
                 0, -minus_sin_theta, cos_theta, 0,
                 0, 0, 0, 1);

    //  绕y轴的旋转矩阵
    cos_theta = Fast_Cos(theta_y);
    minus_sin_theta = -Fast_Sin(theta_y);
    Mat_Init_4X4(&my_inv,
                 cos_theta, 0, -minus_sin_theta, 0,
                 0, 1, 0, 0,
                 minus_sin_theta, 0, cos_theta, 0,
                 0, 0, 0, 1);

    //  绕z轴的旋转矩阵
    cos_theta = Fast_Cos(theta_z);
    minus_sin_theta = -Fast_Sin(theta_z);
    Mat_Init_4X4(&mz_inv,
                 cos_theta, minus_sin_theta, 0, 0,
                 -minus_sin_theta, cos_theta, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1);

    // 根据不同的顺序来进行三个旋转矩阵的乘法
    switch (cam_rot_seq)
    {
    case CAM_ROT_SEQ_XYZ:
    {
        Mat_Mul_4X4(&mx_inv, &my_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mz_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_YXZ:
    {
        Mat_Mul_4X4(&my_inv, &mx_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mz_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_XZY:
    {
        Mat_Mul_4X4(&mx_inv, &mz_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &my_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_YZX:
    {
        Mat_Mul_4X4(&my_inv, &mz_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mx_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_ZYX:
    {
        Mat_Mul_4X4(&mz_inv, &my_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &mx_inv, &mrot);
    }
    break;

    case CAM_ROT_SEQ_ZXY:
    {
        Mat_Mul_4X4(&mz_inv, &mx_inv, &mtmp);
        Mat_Mul_4X4(&mtmp, &my_inv, &mrot);
    }
    break;

    default:
        break;
    }

    //将平移矩阵和旋转矩阵相乘的结果，保存到cam->mcam
    Mat_Mul_4X4(&mt_inv, &mrot, &cam->mcam);
}

////////////////////////////////////////////////////////////

//加载plg格式文件
int Load_OBJECT4DV1_PLG(OBJECT4DV1_PTR obj, // pointer to object
                        char *filename,     // filename of plg file
                        VECTOR4D_PTR scale, // initial scaling factors
                        VECTOR4D_PTR pos,   // initial position
                        VECTOR4D_PTR rot)   // initial rotations
{
    // this function loads a plg object in off disk, additionally
    // it allows the caller to scale, position, and rotate the object
    // to save extra calls later for non-dynamic objects

    FILE *fp;         // file pointer
    char buffer[256]; // working buffer

    char *token_string; // pointer to actual token text, ready for parsing

    // file format review, note types at end of each description
    // # this is a comment

    // # object descriptor
    // object_name_string num_verts_int num_polys_int

    // # vertex list
    // x0_float y0_float z0_float
    // x1_float y1_float z1_float
    // x2_float y2_float z2_float
    // .
    // .
    // xn_float yn_float zn_float
    //
    // # polygon list
    // surface_description_ushort num_verts_int v0_index_int v1_index_int ..  vn_index_int
    // .
    // .
    // surface_description_ushort num_verts_int v0_index_int v1_index_int ..  vn_index_int

    // lets keep it simple and assume one element per line
    // hence we have to find the object descriptor, read it in, then the
    // vertex list and read it in, and finally the polygon list -- simple :)

    // Step 1: clear out the object and initialize it a bit
    memset(obj, 0, sizeof(OBJECT4DV1));

    // set state of object to active and visible
    obj->state = OBJECT4DV1_STATE_ACTIVE | OBJECT4DV1_STATE_VISIBLE;

    // set position of object
    obj->world_pos.x = pos->x;
    obj->world_pos.y = pos->y;
    obj->world_pos.z = pos->z;
    obj->world_pos.w = pos->w;

    // Step 2: open the file for reading
    if (!(fp = fopen(filename, "r")))
    {
        // Write_Error("Couldn't open PLG file %s.", filename);
        return (0);
    } // end if

    // Step 3: get the first token string which should be the object descriptor
    if (!(token_string = Get_Line_PLG(buffer, 255, fp)))
    {
        // Write_Error("PLG file error with file %s (object descriptor invalid).", filename);
        return (0);
    } // end if

    // Write_Error("Object Descriptor: %s", token_string);

    // parse out the info object
    sscanf(token_string, "%s %d %d", obj->name, &obj->num_vertices, &obj->num_polys);

    // Step 4: load the vertex list
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // get the next vertex
        if (!(token_string = Get_Line_PLG(buffer, 255, fp)))
        {
            // Write_Error("PLG file error with file %s (vertex list invalid).", filename);
            return (0);
        } // end if

        // parse out vertex
        sscanf(token_string, "%f %f %f", &obj->vlist_local[vertex].x,
               &obj->vlist_local[vertex].y,
               &obj->vlist_local[vertex].z);
        obj->vlist_local[vertex].w = 1;

        // scale vertices
        obj->vlist_local[vertex].x *= scale->x;
        obj->vlist_local[vertex].y *= scale->y;
        obj->vlist_local[vertex].z *= scale->z;

        // Write_Error("\nVertex %d = %f, %f, %f, %f", vertex,
        //             obj->vlist_local[vertex].x,
        //             obj->vlist_local[vertex].y,
        //             obj->vlist_local[vertex].z,
        //             obj->vlist_local[vertex].w);

    } // end for vertex

    // compute average and max radius
    Compute_OBJECT4DV1_Radius(obj);

    // Write_Error("\nObject average radius = %f, max radius = %f",
    //             obj->avg_radius, obj->max_radius);

    int poly_surface_desc = 0; // PLG/PLX surface descriptor
    int poly_num_verts = 0;    // number of vertices for current poly (always 3)
    char tmp_string[8];        // temp string to hold surface descriptor in and
                               // test if it need to be converted from hex

    // Step 5: load the polygon list
    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        // get the next polygon descriptor
        if (!(token_string = Get_Line_PLG(buffer, 255, fp)))
        {
            // Write_Error("PLG file error with file %s (polygon descriptor invalid).", filename);
            return (0);
        } // end if

        // Write_Error("\nPolygon %d:", poly);

        // each vertex list MUST have 3 vertices since we made this a rule that all models
        // must be constructed of triangles
        // read in surface descriptor, number of vertices, and vertex list
        sscanf(token_string, "%s %d %d %d %d", tmp_string,
               &poly_num_verts, // should always be 3
               &obj->plist[poly].vert[0],
               &obj->plist[poly].vert[1],
               &obj->plist[poly].vert[2]);

        // since we are allowing the surface descriptor to be in hex format
        // with a leading "0x" we need to test for it
        if (tmp_string[0] == '0' && toupper(tmp_string[1]) == 'X')
            sscanf(tmp_string, "%x", &poly_surface_desc);
        else
            poly_surface_desc = atoi(tmp_string);

        // point polygon vertex list to object's vertex list
        // note that this is redundant since the polylist is contained
        // within the object in this case and its up to the user to select
        // whether the local or transformed vertex list is used when building up
        // polygon geometry, might be a better idea to set to NULL in the context
        // of polygons that are part of an object
        obj->plist[poly].vlist = obj->vlist_local;

        // Write_Error("\nSurface Desc = 0x%.4x, num_verts = %d, vert_indices [%d, %d, %d]",
        //             poly_surface_desc,
        //             poly_num_verts,
        //             obj->plist[poly].vert[0],
        //             obj->plist[poly].vert[1],
        //             obj->plist[poly].vert[2]);

        // now we that we have the vertex list and we have entered the polygon
        // vertex index data into the polygon itself, now let's analyze the surface
        // descriptor and set the fields for the polygon based on the description

        // extract out each field of data from the surface descriptor
        // first let's get the single/double sided stuff out of the way
        if ((poly_surface_desc & PLX_2SIDED_FLAG))
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_2SIDED);
            // Write_Error("\n2 sided.");
        } // end if
        else
        {
            // one sided
            // Write_Error("\n1 sided.");
        } // end else

        // now let's set the color type and color
        if ((poly_surface_desc & PLX_COLOR_MODE_RGB_FLAG))
        {
            // this is an RGB 4.4.4 surface
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_RGB16);

            // now extract color and copy into polygon color
            // field in proper 16-bit format
            // 0x0RGB is the format, 4 bits per pixel
            int red = ((poly_surface_desc & 0x0f00) >> 8);
            int green = ((poly_surface_desc & 0x00f0) >> 4);
            int blue = (poly_surface_desc & 0x000f);

            // although the data is always in 4.4.4 format, the graphics card
            // is either 5.5.5 or 5.6.5, but our virtual color system translates
            // 8.8.8 into 5.5.5 or 5.6.5 for us, but we have to first scale all
            // these 4.4.4 values into 8.8.8
            //obj->plist[poly].color = RGB16Bit(red * 16, green * 16, blue * 16);

            //obj->plist[poly].color = RGB16Bit565(red * 16, green * 16, blue * 16);
            obj->plist[poly].color = RGB16Bit565(red * 16, green * 16, blue * 16);

            // IUINT32 c = (red << 16) | (green << 8) | blue;
            // obj->plist[poly].color = c;

            // Write_Error("\nRGB color = [%d, %d, %d]", red, green, blue);
        } // end if
        else
        {
            // this is an 8-bit color indexed surface
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_8BITCOLOR);

            // and simple extract the last 8 bits and that's the color index
            obj->plist[poly].color = (poly_surface_desc & 0x00ff);

            // Write_Error("\n8-bit color index = %d", obj->plist[poly].color);

        } // end else

        // handle shading mode
        int shade_mode = (poly_surface_desc & PLX_SHADE_MODE_MASK);

        // set polygon shading mode
        switch (shade_mode)
        {
        case PLX_SHADE_MODE_PURE_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_PURE);
            // Write_Error("\nShade mode = pure");
        }
        break;

        case PLX_SHADE_MODE_FLAT_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_FLAT);
            // Write_Error("\nShade mode = flat");
        }
        break;

        case PLX_SHADE_MODE_GOURAUD_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_GOURAUD);
            // Write_Error("\nShade mode = gouraud");
        }
        break;

        case PLX_SHADE_MODE_PHONG_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_PHONG);
            // Write_Error("\nShade mode = phong");
        }
        break;

        default:
            break;
        } // end switch

        // finally set the polygon to active
        obj->plist[poly].state = POLY4DV1_STATE_ACTIVE;

    } // end for poly

    // close the file
    fclose(fp);

    // return success
    return (1);

} // end Load_OBJECT4DV1_PLG

char *Get_Line_PLG(char *buffer, int maxlength, FILE *fp)
{
    // this little helper function simply read past comments
    // and blank lines in a PLG file and always returns full
    // lines with something on them on NULL if the file is empty

    int index = 0;  // general index
    int length = 0; // general length

    // enter into parsing loop
    while (1)
    {
        // read the next line
        if (!fgets(buffer, maxlength, fp))
            return (NULL);

        // kill the whitespace
        for (length = strlen(buffer), index = 0; isspace(buffer[index]); index++)
            ;

        // test if this was a blank line or a comment
        if (index >= length || buffer[index] == '#')
            continue;

        // at this point we have a good line
        return (&buffer[index]);
    } // end while

} // end Get_Line_PLG

float Compute_OBJECT4DV1_Radius(OBJECT4DV1_PTR obj)
{
    // this function computes the average and maximum radius for
    // sent object and opdates the object data

    // reset incase there's any residue
    obj->avg_radius = 0;
    obj->max_radius = 0;

    // loop thru and compute radius
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // update the average and maximum radius
        float dist_to_vertex =
            sqrt(obj->vlist_local[vertex].x * obj->vlist_local[vertex].x +
                 obj->vlist_local[vertex].y * obj->vlist_local[vertex].y +
                 obj->vlist_local[vertex].z * obj->vlist_local[vertex].z);

        // accumulate total radius
        obj->avg_radius += dist_to_vertex;

        // update maximum radius
        if (dist_to_vertex > obj->max_radius)
            obj->max_radius = dist_to_vertex;

    } // end for vertex

    // finallize average radius computation
    obj->avg_radius /= obj->num_vertices;

    // return max radius
    return (obj->max_radius);

} // end Compute_OBJECT4DV1_Radius

void VECTOR4D_Build(VECTOR4D_PTR init, VECTOR4D_PTR term, VECTOR4D_PTR result)
{
    // build a 4d vector
    result->x = term->x - init->x;
    result->y = term->y - init->y;
    result->z = term->z - init->z;
    result->w = 1;

} // end VECTOR4D_Build

//向量点乘
float VECTOR4D_Dot(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    return ((va->x * vb->x) + (va->y * vb->y) + (va->z * vb->z));
}

//向量叉乘
void VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vn)
{
    vn->x = ((va->y * vb->z) - (va->z * vb->y));
    vn->y = -((va->x * vb->z) - (va->z * vb->x));
    vn->z = ((va->x * vb->y) - (va->y * vb->x));
    vn->w = 1;
}

//向量叉乘
VECTOR4D VECTOR4D_Cross(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    VECTOR4D vn;
    vn.x = ((va->y * vb->z) - (va->z * vb->y));
    vn.y = -((va->x * vb->z) - (va->z * vb->x));
    vn.z = ((va->x * vb->y) - (va->y * vb->x));
    vn.w = 1;
    return (vn);
}

int Light_RENDERLIST4DV1_World16(RENDERLIST4DV1_PTR rend_list, // list to process
                                 CAM4DV1_PTR cam,              // camera position
                                 LIGHTV1_PTR lights,           // light list (might have more than one)
                                 int max_lights)               // maximum lights in list
{
    // 16-bit version of function
    // function lights the enture rendering list based on the sent lights and camera. the function supports
    // constant/pure shading (emmisive), flat shading with ambient, infinite, point lights, and spot lights
    // note that this lighting function is rather brute force and simply follows the math, however
    // there are some clever integer operations that are used in scale 256 rather than going to floating
    // point, but why? floating point and ints are the same speed, HOWEVER, the conversion to and from floating
    // point can be cycle intensive, so if you can keep your calcs in ints then you can gain some speed
    // also note, type 1 spot lights are simply point lights with direction, the "cone" is more of a function
    // of the falloff due to attenuation, but they still look like spot lights
    // type 2 spot lights are implemented with the intensity having a dot product relationship with the
    // angle from the surface point to the light direction just like in the optimized model, but the pf term
    // that is used for a concentration control must be 1,2,3,.... integral and non-fractional

    // also note since we are dealing with a rendering list and not object, the final lit color is
    // immediately written over the real color

    unsigned int r_base, g_base, b_base, // base color being lit
        r_sum, g_sum, b_sum,             // sum of lighting process over all lights
        shaded_color;                    // final color

    float dp,  // dot product
        dist,  // distance from light to surface
        i,     // general intensities
        nl,    // length of normal
        atten; // attenuation computations

    // for each valid poly, light it...
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire polygon
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // light this polygon if and only if it's not clipped, not culled,
        // active, and visible
        if (!(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue; // move onto next poly

        // we will use the transformed polygon vertex list since the backface removal
        // only makes sense at the world coord stage further of the pipeline

        // test the lighting mode of the polygon (use flat for flat, gouraud))
        if (curr_poly->attr & POLY4DV1_ATTR_SHADE_MODE_FLAT || curr_poly->attr & POLY4DV1_ATTR_SHADE_MODE_GOURAUD)
        {
            // step 1: extract the base color out in RGB mode
            //if (dd_pixel_format == DD_PIXEL_FORMAT565)
            if (1)
            {
                _RGB565FROM16BIT(curr_poly->color, &r_base, &g_base, &b_base);

                // scale to 8 bit
                r_base <<= 3;
                g_base <<= 2;
                b_base <<= 3;
            } // end if
            else
            {
                // _RGB555FROM16BIT(curr_poly->color, &r_base, &g_base, &b_base);

                // // scale to 8 bit
                // r_base <<= 3;
                // g_base <<= 3;
                // b_base <<= 3;
            } // end if

            // initialize color sum
            r_sum = 0;
            g_sum = 0;
            b_sum = 0;

            // loop thru lights
            for (int curr_light = 0; curr_light < max_lights; curr_light++)
            {
                // is this light active
                if (lights[curr_light].state == LIGHTV1_STATE_OFF)
                    continue;

                // what kind of light are we dealing with
                if (lights[curr_light].attr & LIGHTV1_ATTR_AMBIENT)
                {
                    // simply multiply each channel against the color of the
                    // polygon then divide by 256 to scale back to 0..255
                    // use a shift in real life!!! >> 8
                    r_sum += ((lights[curr_light].c_ambient.r * r_base) / 256);
                    g_sum += ((lights[curr_light].c_ambient.g * g_base) / 256);
                    b_sum += ((lights[curr_light].c_ambient.b * b_base) / 256);

                    // there better only be one ambient light!

                } // end if
                else if (lights[curr_light].attr & LIGHTV1_ATTR_INFINITE)
                {
                    // infinite lighting, we need the surface normal, and the direction
                    // of the light source

                    // we need to compute the normal of this polygon face, and recall
                    // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
                    VECTOR4D u, v, n;

                    // build u, v
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[1], &u);
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[2], &v);

                    // compute cross product
                    VECTOR4D_Cross(&u, &v, &n);

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    nl = VECTOR4D_Length_Fast(&n);

                    // ok, recalling the lighting model for infinite lights
                    // I(d)dir = I0dir * Cldir
                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    dp = VECTOR4D_Dot(&n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp / nl;
                        r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                } // end if infinite light
                else if (lights[curr_light].attr & LIGHTV1_ATTR_POINT)
                {
                    // perform point light computations
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    // we need to compute the normal of this polygon face, and recall
                    // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
                    VECTOR4D u, v, n, l;

                    // build u, v
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[1], &u);
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[2], &v);

                    // compute cross product
                    VECTOR4D_Cross(&u, &v, &n);

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    nl = VECTOR4D_Length_Fast(&n);

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0], &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles
                    dp = VECTOR4D_Dot(&n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (nl * dist * atten);

                        r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                } // end if point
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT1)
                {
                    // perform spotlight/point computations simplified model that uses
                    // point light WITH a direction to simulate a spotlight
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    // we need to compute the normal of this polygon face, and recall
                    // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
                    VECTOR4D u, v, n, l;

                    // build u, v
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[1], &u);
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[2], &v);

                    // compute cross product (we need -n, so do vxu)
                    VECTOR4D_Cross(&v, &u, &n);

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    nl = VECTOR4D_Length_Fast(&n);

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0], &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // note that I use the direction of the light here rather than a the vector to the light
                    // thus we are taking orientation into account which is similar to the spotlight model
                    dp = VECTOR4D_Dot(&n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (nl * atten);

                        r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                }                                                           // end if spotlight1
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT2) // simple version
                {
                    // perform spot light computations
                    // light model for spot light simple version is once again:
                    //         	     I0spotlight * Clspotlight * MAX( (l . s), 0)^pf
                    // I(d)spotlight = __________________________________________
                    //               		 kc + kl*d + kq*d2
                    // Where d = |p - s|, and pf = power factor

                    // thus it's almost identical to the point, but has the extra term in the numerator
                    // relating the angle between the light source and the point on the surface

                    // we need to compute the normal of this polygon face, and recall
                    // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
                    VECTOR4D u, v, n, d, s;

                    // build u, v
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[1], &u);
                    VECTOR4D_Build(&curr_poly->tvlist[0], &curr_poly->tvlist[2], &v);

                    // compute cross product (v x u, to invert n)
                    VECTOR4D_Cross(&v, &u, &n);

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    nl = VECTOR4D_Length_Fast(&n);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles
                    dp = VECTOR4D_Dot(&n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[0], &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dist = VECTOR4D_Length_Fast(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dist;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (nl * atten);

                            r_sum += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);

                        } // end if

                    } // end if

                } // end if spot light

            } // end for light

            // make sure colors aren't out of range
            if (r_sum > 255)
                r_sum = 255;
            if (g_sum > 255)
                g_sum = 255;
            if (b_sum > 255)
                b_sum = 255;

            // write the color over current color
            // IUINT32 c = (r_sum << 16) | (g_sum << 8) | b_sum;
            // curr_poly->color = c;

            curr_poly->color = RGB16Bit565(r_sum, g_sum, b_sum);
            //curr_poly->color = RGB16Bit(r_sum, g_sum, b_sum);

        }    // end if
        else // assume POLY4DV1_ATTR_SHADE_MODE_CONSTANT
        {
            // emmisive shading only, do nothing
            // ...
        } // end if

    } // end for poly

    // return success
    return (1);

} // end Light_RENDERLIST4DV1_World16

int num_lights;             // current number of lights
LIGHTV1 lights[MAX_LIGHTS]; // lights in system

int Reset_Lights_LIGHTV1(void)
{
    // this function simply resets all lights in the system
    static int first_time = 1;

    memset(lights, 0, MAX_LIGHTS * sizeof(LIGHTV1));

    // reset number of lights
    num_lights = 0;

    // reset first time
    first_time = 0;

    // return success
    return (1);

} // end Reset_Lights_LIGHTV1

int Init_Light_LIGHTV1(int index,          // index of light to create (0..MAX_LIGHTS-1)
                       int _state,         // state of light
                       int _attr,          // type of light, and extra qualifiers
                       RGBAV1 _c_ambient,  // ambient light intensity
                       RGBAV1 _c_diffuse,  // diffuse light intensity
                       RGBAV1 _c_specular, // specular light intensity
                       POINT4D_PTR _pos,   // position of light
                       VECTOR4D_PTR _dir,  // direction of light
                       float _kc,          // attenuation factors
                       float _kl,
                       float _kq,
                       float _spot_inner, // inner angle for spot light
                       float _spot_outer, // outer angle for spot light
                       float _pf)         // power factor/falloff for spot lights
{
    // this function initializes a light based on the flags sent in _attr, values that
    // aren't needed are set to 0 by caller

    // make sure light is in range
    if (index < 0 || index >= MAX_LIGHTS)
        return (0);

    // all good, initialize the light (many fields may be dead)
    lights[index].state = _state; // state of light
    lights[index].id = index;     // id of light
    lights[index].attr = _attr;   // type of light, and extra qualifiers

    lights[index].c_ambient = _c_ambient;   // ambient light intensity
    lights[index].c_diffuse = _c_diffuse;   // diffuse light intensity
    lights[index].c_specular = _c_specular; // specular light intensity

    lights[index].kc = _kc; // constant, linear, and quadratic attenuation factors
    lights[index].kl = _kl;
    lights[index].kq = _kq;

    if (_pos)
        VECTOR4D_COPY(&lights[index].pos, _pos); // position of light

    if (_dir)
    {
        VECTOR4D_COPY(&lights[index].dir, _dir); // direction of light
        // normalize it
        VECTOR4D_Normalize(&lights[index].dir);

    } // end if

    lights[index].spot_inner = _spot_inner; // inner angle for spot light
    lights[index].spot_outer = _spot_outer; // outer angle for spot light
    lights[index].pf = _pf;                 // power factor/falloff for spot lights

    // return light index as success
    return (index);

} // end Create_Light_LIGHTV1

void VECTOR4D_Normalize(VECTOR4D_PTR va)
{
    // normalizes the sent vector and returns the result

    // compute length
    float length = sqrtf(va->x * va->x + va->y * va->y + va->z * va->z);

    // test for zero length vector
    // if found return zero vector
    if (length < EPSILON_E5)
        return;

    float length_inv = 1.0 / length;

    // compute normalized version of vector
    va->x *= length_inv;
    va->y *= length_inv;
    va->z *= length_inv;
    va->w = 1;

} // end VECTOR4D_Normalize

float VECTOR4D_Length_Fast(VECTOR4D_PTR va)
{
    // computes the magnitude of a vector using an approximation
    // very fast
    return (Fast_Distance_3D(va->x, va->y, va->z));

} // end VECTOR4D_Length_Fast

float Fast_Distance_3D(float fx, float fy, float fz)
{
    // this function computes the distance from the origin to x,y,z

    int temp;    // used for swaping
    int x, y, z; // used for algorithm

    // make sure values are all positive
    x = fabs(fx) * 1024;
    y = fabs(fy) * 1024;
    z = fabs(fz) * 1024;

    // sort values
    if (y < x)
        SWAP(x, y, temp)

    if (z < y)
        SWAP(y, z, temp)

    if (y < x)
        SWAP(x, y, temp)

    int dist = (z + 11 * (y >> 5) + (x >> 2));

    // compute distance with 8% error
    return ((float)(dist >> 10));

} // end Fast_Distance_3D

void Sort_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, int sort_method)
{
    // this function sorts the rendering list based on the polygon z-values
    // the specific sorting method is controlled by sending in control flags
    // #define SORT_POLYLIST_AVGZ  0 - sorts on average of all vertices
    // #define SORT_POLYLIST_NEARZ 1 - sorts on closest z vertex of each poly
    // #define SORT_POLYLIST_FARZ  2 - sorts on farthest z vertex of each poly

    switch (sort_method)
    {
    case SORT_POLYLIST_AVGZ: //  - sorts on average of all vertices
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV1_PTR), Compare_AvgZ_POLYF4DV1);
    }
    break;

    case SORT_POLYLIST_NEARZ: // - sorts on closest z vertex of each poly
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV1_PTR), Compare_NearZ_POLYF4DV1);
    }
    break;

    case SORT_POLYLIST_FARZ: //  - sorts on farthest z vertex of each poly
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV1_PTR), Compare_FarZ_POLYF4DV1);
    }
    break;

    default:
        break;
    } // end switch

} // end Sort_RENDERLIST4DV1

int Compare_AvgZ_POLYF4DV1(const void *arg1, const void *arg2)
{
    // this function comapares the average z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV1_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV1_PTR *)(arg1));
    poly_2 = *((POLYF4DV1_PTR *)(arg2));

    // compute z average of each polygon
    z1 = (float)0.33333 * (poly_1->tvlist[0].z + poly_1->tvlist[1].z + poly_1->tvlist[2].z);

    // now polygon 2
    z2 = (float)0.33333 * (poly_2->tvlist[0].z + poly_2->tvlist[1].z + poly_2->tvlist[2].z);

    // compare z1 and z2, such that polys' will be sorted in descending Z order
    if (z1 > z2)
        return (-1);
    else if (z1 < z2)
        return (1);
    else
        return (0);

} // end Compare_AvgZ_POLYF4DV1

int Compare_NearZ_POLYF4DV1(const void *arg1, const void *arg2)
{
    // this function comapares the closest z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV1_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV1_PTR *)(arg1));
    poly_2 = *((POLYF4DV1_PTR *)(arg2));

    // compute the near z of each polygon
    z1 = MIN(poly_1->tvlist[0].z, poly_1->tvlist[1].z);
    z1 = MIN(z1, poly_1->tvlist[2].z);

    z2 = MIN(poly_2->tvlist[0].z, poly_2->tvlist[1].z);
    z2 = MIN(z2, poly_2->tvlist[2].z);

    // compare z1 and z2, such that polys' will be sorted in descending Z order
    if (z1 > z2)
        return (-1);
    else if (z1 < z2)
        return (1);
    else
        return (0);

} // end Compare_NearZ_POLYF4DV1

////////////////////////////////////////////////////////////////////////////////

int Compare_FarZ_POLYF4DV1(const void *arg1, const void *arg2)
{
    // this function comapares the farthest z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV1_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV1_PTR *)(arg1));
    poly_2 = *((POLYF4DV1_PTR *)(arg2));

    // compute the near z of each polygon
    z1 = MAX(poly_1->tvlist[0].z, poly_1->tvlist[1].z);
    z1 = MAX(z1, poly_1->tvlist[2].z);

    z2 = MAX(poly_2->tvlist[0].z, poly_2->tvlist[1].z);
    z2 = MAX(z2, poly_2->tvlist[2].z);

    // compare z1 and z2, such that polys' will be sorted in descending Z order
    if (z1 > z2)
        return (-1);
    else if (z1 < z2)
        return (1);
    else
        return (0);

} // end Compare_FarZ_POLYF4DV1