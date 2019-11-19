#include "datastructure.h"

//extern 变量
int num_lights;             // current number of lights
LIGHTV1 lights[MAX_LIGHTS]; // lights in system

MATV1 materials[MAX_MATERIALS]; // materials in system
int num_materials;              // current number of materials

char texture_path[80] = "./"; // root path to ALL textures, make current directory for now

BITMAP_FILE bitmap16bit; // a 16 bit bitmap file
// these are overwritten globally by DDraw_Init()

int screen_width = SCREEN_WIDTH,   // width of screen
    screen_height = SCREEN_HEIGHT, // height of screen
    screen_bpp = SCREEN_BPP,       // bits per pixel
    screen_windowed = 0;           // is this a windowed app?

PALETTEENTRY palette[MAX_COLORS_PALETTE]; // color palette

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
        // //Write_Error("Couldn't open PLG file %s.", filename);
        return (0);
    } // end if

    // Step 3: get the first token string which should be the object descriptor
    if (!(token_string = Get_Line_PLG(buffer, 255, fp)))
    {
        // //Write_Error("PLG file error with file %s (object descriptor invalid).", filename);
        return (0);
    } // end if

    // //Write_Error("Object Descriptor: %s", token_string);

    // parse out the info object
    sscanf(token_string, "%s %d %d", obj->name, &obj->num_vertices, &obj->num_polys);

    // Step 4: load the vertex list
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // get the next vertex
        if (!(token_string = Get_Line_PLG(buffer, 255, fp)))
        {
            // //Write_Error("PLG file error with file %s (vertex list invalid).", filename);
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

        // //Write_Error("\nVertex %d = %f, %f, %f, %f", vertex,
        //             obj->vlist_local[vertex].x,
        //             obj->vlist_local[vertex].y,
        //             obj->vlist_local[vertex].z,
        //             obj->vlist_local[vertex].w);

    } // end for vertex

    // compute average and max radius
    Compute_OBJECT4DV1_Radius(obj);

    // //Write_Error("\nObject average radius = %f, max radius = %f",
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
            // //Write_Error("PLG file error with file %s (polygon descriptor invalid).", filename);
            return (0);
        } // end if

        // //Write_Error("\nPolygon %d:", poly);

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

        // //Write_Error("\nSurface Desc = 0x%.4x, num_verts = %d, vert_indices [%d, %d, %d]",
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
            // //Write_Error("\n2 sided.");
        } // end if
        else
        {
            // one sided
            // //Write_Error("\n1 sided.");
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

            // //Write_Error("\nRGB color = [%d, %d, %d]", red, green, blue);
        } // end if
        else
        {
            // this is an 8-bit color indexed surface
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_8BITCOLOR);

            // and simple extract the last 8 bits and that's the color index
            obj->plist[poly].color = (poly_surface_desc & 0x00ff);

            // //Write_Error("\n8-bit color index = %d", obj->plist[poly].color);

        } // end else

        // handle shading mode
        int shade_mode = (poly_surface_desc & PLX_SHADE_MODE_MASK);

        // set polygon shading mode
        switch (shade_mode)
        {
        case PLX_SHADE_MODE_PURE_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_PURE);
            // //Write_Error("\nShade mode = pure");
        }
        break;

        case PLX_SHADE_MODE_FLAT_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_FLAT);
            // //Write_Error("\nShade mode = flat");
        }
        break;

        case PLX_SHADE_MODE_GOURAUD_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_GOURAUD);
            // //Write_Error("\nShade mode = gouraud");
        }
        break;

        case PLX_SHADE_MODE_PHONG_FLAG:
        {
            SET_BIT(obj->plist[poly].attr, POLY4DV1_ATTR_SHADE_MODE_PHONG);
            // //Write_Error("\nShade mode = phong");
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

int Load_OBJECT4DV2_COB(OBJECT4DV2_PTR obj, // pointer to object
                        char *filename,     // filename of Caligari COB file
                        VECTOR4D_PTR scale, // initial scaling factors
                        VECTOR4D_PTR pos,   // initial position
                        VECTOR4D_PTR rot,   // initial rotations
                        int vertex_flags)   // flags to re-order vertices
                                            // and perform transforms
{
    // this function loads a Caligari TrueSpace .COB file object in off disk, additionally
    // it allows the caller to scale, position, and rotate the object
    // to save extra calls later for non-dynamic objects, note that this function
    // works with a OBJECT4DV2 which has support for textures, but not materials, etc,
    // however we will still parse out the material stuff and get them ready for the
    // next incarnation objects, so we can re-use this code to support those features
    // also, since this version IS going to read in the texture map and texture coordinates
    // we have a couple issues to think about, first COB format like absolute texture paths
    // we can't have that, so we will simple extract out ONLY the texture map bitmap name
    // and use the global texture path variable to build a real file path, also texture
    // coordinates are in 0..1 0..1 form, I still haven't decided if I want to use absolute
    // coordinates or 0..1 0..1, but right now the affine texture mapper uses

    // create a parser object
    CPARSERV1 parser;

    char seps[16];          // seperators for token scanning
    char token_buffer[256]; // used as working buffer for token
    char *token;            // pointer to next token

    int r, g, b; // working colors

    // cache for texture vertices
    VERTEX2DF texture_vertices[OBJECT4DV2_MAX_VERTICES];

    int num_texture_vertices = 0;

    MATRIX4X4 mat_local, // storage for local transform if user requests it in cob format
        mat_world;       // "   " for local to world " "

    // initialize matrices
    MAT_IDENTITY_4X4(&mat_local);
    MAT_IDENTITY_4X4(&mat_world);

    // Step 1: clear out the object and initialize it a bit
    memset(obj, 0, sizeof(OBJECT4DV2));

    // set state of object to active and visible
    obj->state = OBJECT4DV2_STATE_ACTIVE | OBJECT4DV2_STATE_VISIBLE;

    // set number of frames
    obj->num_frames = 1;
    obj->curr_frame = 0;
    obj->attr = OBJECT4DV2_ATTR_SINGLE_FRAME;

    // set position of object is caller requested position
    if (pos)
    {
        // set position of object
        obj->world_pos.x = pos->x;
        obj->world_pos.y = pos->y;
        obj->world_pos.z = pos->z;
        obj->world_pos.w = pos->w;
    } // end
    else
    {
        // set it to (0,0,0,1)
        obj->world_pos.x = 0;
        obj->world_pos.y = 0;
        obj->world_pos.z = 0;
        obj->world_pos.w = 1;
    } // end else

    // Step 2: open the file for reading using the parser
    if (!parser.Open(filename))
    {
        //Write_Error("Couldn't open .COB file %s.", filename);
        return (0);
    } // end if

    // Step 3:

    // lets find the name of the object first
    while (1)
    {
        // get the next line, we are looking for "Name"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("Image 'name' not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Name'] [s>0]"))
        {
            // name should be in second string variable, index 1
            strcpy(obj->name, parser.pstrings[1]);
            //Write_Error("\nCOB Reader Object Name: %s", obj->name);

            break;
        } // end if

    } // end while

    // step 4: get local and world transforms and store them

    // center 0 0 0
    // x axis 1 0 0
    // y axis 0 1 0
    // z axis 0 0 1

    while (1)
    {
        // get the next line, we are looking for "center"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("Center not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['center'] [f] [f] [f]"))
        {
            // the "center" holds the translation factors, so place in
            // last row of homogeneous matrix, note that these are row vectors
            // that we need to drop in each column of matrix
            mat_local.M[3][0] = -parser.pfloats[0]; // center x
            mat_local.M[3][1] = -parser.pfloats[1]; // center y
            mat_local.M[3][2] = -parser.pfloats[2]; // center z

            // ok now, the next 3 lines should be the x,y,z transform vectors
            // so build up

            // "x axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "['x'] ['axis'] [f] [f] [f]");

            // place row in x column of transform matrix
            mat_local.M[0][0] = parser.pfloats[0]; // rxx
            mat_local.M[1][0] = parser.pfloats[1]; // rxy
            mat_local.M[2][0] = parser.pfloats[2]; // rxz

            // "y axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "['y'] ['axis'] [f] [f] [f]");

            // place row in y column of transform matrix
            mat_local.M[0][1] = parser.pfloats[0]; // ryx
            mat_local.M[1][1] = parser.pfloats[1]; // ryy
            mat_local.M[2][1] = parser.pfloats[2]; // ryz

            // "z axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "['z'] ['axis'] [f] [f] [f]");

            // place row in z column of transform matrix
            mat_local.M[0][2] = parser.pfloats[0]; // rzx
            mat_local.M[1][2] = parser.pfloats[1]; // rzy
            mat_local.M[2][2] = parser.pfloats[2]; // rzz

            Print_Mat_4X4(&mat_local, "Local COB Matrix:");

            break;
        } // end if

    } // end while

    // now "Transform"
    while (1)
    {
        // get the next line, we are looking for "Transform"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("Transform not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Transform']"))
        {

            // "x axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "[f] [f] [f]");

            // place row in x column of transform matrix
            mat_world.M[0][0] = parser.pfloats[0]; // rxx
            mat_world.M[1][0] = parser.pfloats[1]; // rxy
            mat_world.M[2][0] = parser.pfloats[2]; // rxz

            // "y axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "[f] [f] [f]");

            // place row in y column of transform matrix
            mat_world.M[0][1] = parser.pfloats[0]; // ryx
            mat_world.M[1][1] = parser.pfloats[1]; // ryy
            mat_world.M[2][1] = parser.pfloats[2]; // ryz

            // "z axis"
            parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS);
            parser.Pattern_Match(parser.buffer, "[f] [f] [f]");

            // place row in z column of transform matrix
            mat_world.M[0][2] = parser.pfloats[0]; // rzx
            mat_world.M[1][2] = parser.pfloats[1]; // rzy
            mat_world.M[2][2] = parser.pfloats[2]; // rzz

            Print_Mat_4X4(&mat_world, "World COB Matrix:");

            // no need to read in last row, since it's always 0,0,0,1 and we don't use it anyway
            break;

        } // end if

    } // end while

    // step 6: get number of vertices and polys in object
    while (1)
    {
        // get the next line, we are looking for "World Vertices"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("'World Vertices' line not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['World'] ['Vertices'] [i]"))
        {
            // simply extract the number of vertices from the pattern matching
            // output arrays
            obj->num_vertices = parser.pints[0];

            //Write_Error("\nCOB Reader Num Vertices: %d", obj->num_vertices);
            break;

        } // end if

    } // end while

    // allocate the memory for the vertices and number of polys (unknown, so use 3*num_vertices)
    // the call parameters are redundant in this case, but who cares
    if (!Init_OBJECT4DV2(obj, // object to allocate
                         obj->num_vertices,
                         obj->num_vertices * 3,
                         obj->num_frames))
    {
        //Write_Error("\nASC file error with file %s (can't allocate memory).",filename);
    } // end if

    // Step 7: load the vertex list
    // now read in vertex list, format:
    // "d.d d.d d.d"
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // hunt for vertex
        while (1)
        {
            // get the next vertex
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nVertex list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "[f] [f] [f]"))
            {
                // at this point we have the x,y,z in the the pfloats array locations 0,1,2
                obj->vlist_local[vertex].x = parser.pfloats[0];
                obj->vlist_local[vertex].y = parser.pfloats[1];
                obj->vlist_local[vertex].z = parser.pfloats[2];
                obj->vlist_local[vertex].w = 1;

                // do vertex swapping right here, allow muliple swaps, why not!
                // defines for vertex re-ordering flags

                //#define VERTEX_FLAGS_INVERT_X   1    // inverts the Z-coordinates
                //#define VERTEX_FLAGS_INVERT_Y   2    // inverts the Z-coordinates
                //#define VERTEX_FLAGS_INVERT_Z   4    // inverts the Z-coordinates
                //#define VERTEX_FLAGS_SWAP_YZ    8    // transforms a RHS model to a LHS model
                //#define VERTEX_FLAGS_SWAP_XZ    16   // ???
                //#define VERTEX_FLAGS_SWAP_XY    32
                //#define VERTEX_FLAGS_INVERT_WINDING_ORDER 64  // invert winding order from cw to ccw or ccw to cc
                //#define VERTEX_FLAGS_TRANSFORM_LOCAL         512   // if file format has local transform then do it!
                //#define VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD  1024  // if file format has local to world then do it!

                VECTOR4D temp_vector; // temp for calculations

                // now apply local and world transformations encoded in COB format
                if (vertex_flags & VERTEX_FLAGS_TRANSFORM_LOCAL)
                {
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, &mat_local, &temp_vector);
                    VECTOR4D_COPY(&obj->vlist_local[vertex].v, &temp_vector);
                } // end if

                if (vertex_flags & VERTEX_FLAGS_TRANSFORM_LOCAL_WORLD)
                {
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, &mat_world, &temp_vector);
                    VECTOR4D_COPY(&obj->vlist_local[vertex].v, &temp_vector);
                } // end if

                float temp_f; // used for swapping

                // invert signs?
                if (vertex_flags & VERTEX_FLAGS_INVERT_X)
                    obj->vlist_local[vertex].x = -obj->vlist_local[vertex].x;

                if (vertex_flags & VERTEX_FLAGS_INVERT_Y)
                    obj->vlist_local[vertex].y = -obj->vlist_local[vertex].y;

                if (vertex_flags & VERTEX_FLAGS_INVERT_Z)
                    obj->vlist_local[vertex].z = -obj->vlist_local[vertex].z;

                // swap any axes?
                if (vertex_flags & VERTEX_FLAGS_SWAP_YZ)
                    SWAP(obj->vlist_local[vertex].y, obj->vlist_local[vertex].z, temp_f);

                if (vertex_flags & VERTEX_FLAGS_SWAP_XZ)
                    SWAP(obj->vlist_local[vertex].x, obj->vlist_local[vertex].z, temp_f);

                if (vertex_flags & VERTEX_FLAGS_SWAP_XY)
                    SWAP(obj->vlist_local[vertex].x, obj->vlist_local[vertex].y, temp_f);

                // scale vertices
                if (scale)
                {
                    obj->vlist_local[vertex].x *= scale->x;
                    obj->vlist_local[vertex].y *= scale->y;
                    obj->vlist_local[vertex].z *= scale->z;
                } // end if

                //Write_Error("\nVertex %d = %f, %f, %f, %f", vertex,obj->vlist_local[vertex].x, obj->vlist_local[vertex].y, obj->vlist_local[vertex].z,obj->vlist_local[vertex].w);

                // set point field in this vertex, we need that at least
                SET_BIT(obj->vlist_local[vertex].attr, VERTEX4DTV1_ATTR_POINT);

                // found vertex, break out of while for next pass
                break;

            } // end if

        } // end while

    } // end for vertex

    // compute average and max radius
    Compute_OBJECT4DV2_Radius(obj);

    //Write_Error("\nObject average radius = %f, max radius = %f",

    // step 8: get number of texture vertices
    while (1)
    {
        // get the next line, we are looking for "Texture Vertices ddd"
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("'Texture Vertices' line not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Texture'] ['Vertices'] [i]"))
        {
            // simply extract the number of texture vertices from the pattern matching
            // output arrays
            num_texture_vertices = parser.pints[0];

            //Write_Error("\nCOB Reader Texture Vertices: %d", num_texture_vertices);
            break;

        } // end if

    } // end while

    // Step 9: load the texture vertex list in format "U V"
    // "d.d d.d"
    for (int tvertex = 0; tvertex < num_texture_vertices; tvertex++)
    {
        // hunt for texture
        while (1)
        {
            // get the next vertex
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nTexture Vertex list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "[f] [f]"))
            {
                // at this point we have the U V in the the pfloats array locations 0,1 for this
                // texture vertex, store in texture coordinate list
                // note texture coords are in 0..1 format, and must be scaled to texture size
                // after we load the texture
                obj->tlist[tvertex].x = parser.pfloats[0];
                obj->tlist[tvertex].y = parser.pfloats[1];

                //Write_Error("\nTexture Vertex %d: U=%f, V=%f", tvertex, obj->tlist[tvertex].x, obj->tlist[tvertex].y );

                // found vertex, break out of while for next pass
                break;

            } // end if

        } // end while

    } // end for

    // when we load in the polygons then we will copy the texture vertices into the polygon
    // vertices assuming that each vertex has a SINGLE texture coordinate, this means that
    // you must NOT use multiple textures on an object! in other words think "skin" this is
    // inline with Quake II md2 format, in 99% of the cases a single object can be textured
    // with a single skin and the texture coordinates can be unique for each vertex and 1:1

    int poly_material[OBJECT4DV2_MAX_POLYS]; // this holds the material index for each polygon
                                             // we need these indices since when reading the file
                                             // we read the polygons BEFORE the materials, so we need
                                             // this data, so we can go back later and extract the material
                                             // that each poly WAS assigned and get the colors out, since
                                             // objects and polygons do not currently support materials

    int material_index_referenced[MAX_MATERIALS]; // used to track if an index has been used yet as a material
                                                  // reference. since we don't know how many materials, we need
                                                  // a way to count them up, but if we have seen a material reference
                                                  // more than once then we don't increment the total number of materials
                                                  // this array is for this

    // clear out reference array
    memset(material_index_referenced, 0, sizeof(material_index_referenced));

    // step 10: load in the polygons
    // poly list starts off with:
    // "Faces ddd:"
    while (1)
    {
        // get next line
        if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
        {
            //Write_Error("\n'Faces' line not found in .COB file %s.", filename);
            return (0);
        } // end if

        // check for pattern?
        if (parser.Pattern_Match(parser.buffer, "['Faces'] [i]"))
        {
            //Write_Error("\nCOB Reader found face list in .COB file %s.", filename);

            // finally set number of polys
            obj->num_polys = parser.pints[0];

            break;
        } // end if
    }     // end while

    // now read each face in format:
    // Face verts nn flags ff mat mm
    // the nn is the number of vertices, always 3
    // the ff is the flags, unused for now, has to do with holes
    // the mm is the material index number

    int poly_surface_desc = 0;    // ASC surface descriptor/material in this case
    int poly_num_verts = 0;       // number of vertices for current poly (always 3)
    int num_materials_object = 0; // number of materials for this object

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        //Write_Error("\nPolygon %d:", poly);
        // hunt until next face is found
        while (1)
        {
            // get the next polygon face
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nface list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "['Face'] ['verts'] [i] ['flags'] [i] ['mat'] [i]"))
            {
                // at this point we have the number of vertices for the polygon, the flags, and it's material index
                // in the integer output array locations 0,1,2

                // store the material index for this polygon for retrieval later, but make sure adjust the
                // the index to take into consideration that the data in parser.pints[2] is 0 based, and we need
                // an index relative to the entire library, so we simply need to add num_materials to offset the
                // index properly, but we will leave this reference zero based for now... and fix up later
                poly_material[poly] = parser.pints[2];

                // update the reference array
                if (material_index_referenced[poly_material[poly]] == 0)
                {
                    // mark as referenced
                    material_index_referenced[poly_material[poly]] = 1;

                    // increment total number of materials for this object
                    num_materials_object++;
                } // end if

                // test if number of vertices is 3
                if (parser.pints[0] != 3)
                {
                    //Write_Error("\nface not a triangle! in .COB file %s.", filename);
                    return (0);
                } // end if

                // now read out the vertex indices and texture indices format:
                // <vindex0, tindex0>  <vindex1, tindex1> <vindex1, tindex1>
                if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                {
                    //Write_Error("\nface list ended abruptly! in .COB file %s.", filename);
                    return (0);
                } // end if

                // lets replace ",<>" with ' ' to make extraction easy
                ReplaceChars(parser.buffer, parser.buffer, ",<>", ' ');
                parser.Pattern_Match(parser.buffer, "[i] [i] [i] [i] [i] [i]");

                // 0,2,4 holds vertex indices
                // 1,3,5 holds texture indices

                // insert polygon, check for winding order invert
                if (vertex_flags & VERTEX_FLAGS_INVERT_WINDING_ORDER)
                {
                    poly_num_verts = 3;
                    obj->plist[poly].vert[0] = parser.pints[4];
                    obj->plist[poly].vert[1] = parser.pints[2];
                    obj->plist[poly].vert[2] = parser.pints[0];

                    // now copy the texture coordinates into the vertices, this
                    // may not be needed if the polygon doesn't have texture mapping
                    // enabled, etc.,

                    // so here's the deal the texture coordinates that
                    // map to vertex 0,1,2 have indices stored in the odd
                    // numbered pints[] locations, so we simply need to copy
                    // the right texture coordinate into the right vertex
                    obj->plist[poly].text[0] = parser.pints[5];
                    obj->plist[poly].text[1] = parser.pints[3];
                    obj->plist[poly].text[2] = parser.pints[1];

                } // end if
                else
                { // leave winding order alone
                    poly_num_verts = 3;
                    obj->plist[poly].vert[0] = parser.pints[0];
                    obj->plist[poly].vert[1] = parser.pints[2];
                    obj->plist[poly].vert[2] = parser.pints[4];

                    // now copy the texture coordinates into the vertices, this
                    // may not be needed if the polygon doesn't have texture mapping
                    // enabled, etc.,

                    // so here's the deal the texture coordinates that
                    // map to vertex 0,1,2 have indices stored in the odd
                    // numbered pints[] locations, so we simply need to copy
                    // the right texture coordinate into the right vertex
                    obj->plist[poly].text[0] = parser.pints[1];
                    obj->plist[poly].text[1] = parser.pints[3];
                    obj->plist[poly].text[2] = parser.pints[5];

                } // end else

                // point polygon vertex list to object's vertex list
                // note that this is redundant since the polylist is contained
                // within the object in this case and its up to the user to select
                // whether the local or transformed vertex list is used when building up
                // polygon geometry, might be a better idea to set to NULL in the context
                // of polygons that are part of an object
                obj->plist[poly].vlist = obj->vlist_local;

                // set texture coordinate list, this is needed
                obj->plist[poly].tlist = obj->tlist;

                // set polygon to active
                obj->plist[poly].state = POLY4DV2_STATE_ACTIVE;

                // found the face, break out of while for another pass
                break;

            } // end if

        } // end while

    } // end for poly

    // now find materials!!! and we are out of here!
    for (int curr_material = 0; curr_material < num_materials_object; curr_material++)
    {
        // hunt for the material header "mat# ddd"
        while (1)
        {
            // get the next polygon material
            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
            {
                //Write_Error("\nmaterial list ended abruptly! in .COB file %s.", filename);
                return (0);
            } // end if

            // check for pattern?
            if (parser.Pattern_Match(parser.buffer, "['mat#'] [i]"))
            {
                // extract the material that is being defined
                int material_index = parser.pints[0];

                // get color of polygon, although it might be irrelevant for a textured surface
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nRGB color ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // replace the , comma's if there are any with spaces
                    ReplaceChars(parser.buffer, parser.buffer, ",", ' ', 1);

                    // look for "rgb float,float,float"
                    if (parser.Pattern_Match(parser.buffer, "['rgb'] [f] [f] [f]"))
                    {
                        // extract data and store color in material libary
                        // pfloats[] 0,1,2,3, has data
                        materials[material_index + num_materials].color.r = (int)(parser.pfloats[0] * 255 + 0.5);
                        materials[material_index + num_materials].color.g = (int)(parser.pfloats[1] * 255 + 0.5);
                        materials[material_index + num_materials].color.b = (int)(parser.pfloats[2] * 255 + 0.5);

                        break; // while looking for rgb
                    }          // end if

                } // end while

                // extract out lighting constants for the heck of it, they are on a line like this:
                // "alpha float ka float ks float exp float ior float"
                // alpha is transparency           0 - 1
                // ka is ambient coefficient       0 - 1
                // ks is specular coefficient      0 - 1
                // exp is highlight power exponent 0 - 1
                // ior is index of refraction (unused)

                // although our engine will have minimal support for these, we might as well get them
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nmaterial properties ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // look for "alpha float ka float ks float exp float ior float"
                    if (parser.Pattern_Match(parser.buffer, "['alpha'] [f] ['ka'] [f] ['ks'] [f] ['exp'] [f]"))
                    {
                        // extract data and store in material libary
                        // pfloats[] 0,1,2,3, has data
                        materials[material_index + num_materials].color.a = (UCHAR)(parser.pfloats[0] * 255 + 0.5);
                        materials[material_index + num_materials].ka = parser.pfloats[1];
                        materials[material_index + num_materials].kd = 1; // hard code for now
                        materials[material_index + num_materials].ks = parser.pfloats[2];
                        materials[material_index + num_materials].power = parser.pfloats[3];

                        // compute material reflectivities in pre-multiplied format to help engine
                        for (int rgb_index = 0; rgb_index < 3; rgb_index++)
                        {
                            // ambient reflectivity
                            materials[material_index + num_materials].ra.rgba_M[rgb_index] =
                                ((UCHAR)(materials[material_index + num_materials].ka *
                                             (float)materials[material_index + num_materials].color.rgba_M[rgb_index] +
                                         0.5));

                            // diffuse reflectivity
                            materials[material_index + num_materials].rd.rgba_M[rgb_index] =
                                ((UCHAR)(materials[material_index + num_materials].kd *
                                             (float)materials[material_index + num_materials].color.rgba_M[rgb_index] +
                                         0.5));

                            // specular reflectivity
                            materials[material_index + num_materials].rs.rgba_M[rgb_index] =
                                ((UCHAR)(materials[material_index + num_materials].ks *
                                             (float)materials[material_index + num_materials].color.rgba_M[rgb_index] +
                                         0.5));

                        } // end for rgb_index

                        break;
                    } // end if

                } // end while

                // now we need to know the shading model, it's a bit tricky, we need to look for the lines
                // "Shader class: color" first, then after this line is:
                // "Shader name: "xxxxxx" (xxxxxx) "
                // where the xxxxx part will be "plain color" and "plain" for colored polys
                // or "texture map" and "caligari texture"  for textures
                // THEN based on that we hunt for "Shader class: reflectance" which is where the type
                // of shading is encoded, we look for the "Shader name: "xxxxxx" (xxxxxx) " again,
                // and based on it's value we map it to our shading system as follows:
                // "constant" -> MATV1_ATTR_SHADE_MODE_CONSTANT
                // "matte"    -> MATV1_ATTR_SHADE_MODE_FLAT
                // "plastic"  -> MATV1_ATTR_SHADE_MODE_GOURAUD
                // "phong"    -> MATV1_ATTR_SHADE_MODE_FASTPHONG
                // and in the case that in the "color" class, we found a "texture map" then the "shading mode" is
                // "texture map" -> MATV1_ATTR_SHADE_MODE_TEXTURE
                // which must be logically or'ed with the other previous modes

                //  look for the "shader class: color"
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader class ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['class:'] ['color']"))
                    {
                        break;
                    } // end if

                } // end while

                // now look for the shader name for this class
                // Shader name: "plain color" or Shader name: "texture map"
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader name ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // replace the " with spaces
                    ReplaceChars(parser.buffer, parser.buffer, "\"", ' ', 1);

                    // is this a "plain color" poly?
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['name:'] ['plain'] ['color']"))
                    {
                        // not much to do this is default, we need to wait for the reflectance type
                        // to tell us the shading mode

                        break;
                    } // end if

                    // is this a "texture map" poly?
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['name:'] ['texture'] ['map']"))
                    {
                        // set the texture mapping flag in material
                        SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_TEXTURE);

                        // almost done, we need the file name of the darn texture map, its in this format:
                        // file name: string "D:\Source\..\models\textures\wall01.bmp"

                        // of course the filename in the quotes will change
                        // so lets hunt until we find it...
                        while (1)
                        {
                            // get the next line
                            if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                            {
                                //Write_Error("\ncouldnt find texture name! in .COB file %s.", filename);
                                return (0);
                            } // end if

                            // replace the " with spaces
                            ReplaceChars(parser.buffer, parser.buffer, "\"", ' ', 1);

                            // is this the file name?
                            if (parser.Pattern_Match(parser.buffer, "['file'] ['name:'] ['string']"))
                            {
                                // and save the FULL filename (useless though since its the path from the
                                // machine that created it, but later we might want some of the info).
                                // filename and path starts at char position 19, 0 indexed
                                memcpy(materials[material_index + num_materials].texture_file, &parser.buffer[18], strlen(parser.buffer) - 18 + 2);

                                // the OBJECT4DV2 is only allowed a single texture, although we are loading in all
                                // the materials, if this is the first texture map, load it, and set a flag disallowing
                                // any more texture loads for the object
                                if (!obj->texture)
                                {
                                    // step 1: allocate memory for bitmap
                                    obj->texture = (BITMAP_IMAGE_PTR)malloc(sizeof(BITMAP_IMAGE));

                                    // load the texture, just use the final file name and the absolute global
                                    // texture path
                                    char filename[80];
                                    char path_filename[80];
                                    // get the filename
                                    Extract_Filename_From_Path(materials[material_index + num_materials].texture_file, filename);

                                    // build the filename with root path
                                    strcpy(path_filename, texture_path);
                                    strcat(path_filename, filename);

                                    // buffer now holds final texture path and file name
                                    // load the bitmap(8/16 bit)
                                    Load_Bitmap_File(&bitmap16bit, path_filename);

                                    // create a proper size and bitdepth bitmap
                                    Create_Bitmap(obj->texture, 0, 0,
                                                  bitmap16bit.bitmapinfoheader.biWidth,
                                                  bitmap16bit.bitmapinfoheader.biHeight,
                                                  bitmap16bit.bitmapinfoheader.biBitCount);

                                    // load the bitmap image (later make this 8/16 bit)
                                    if (obj->texture->bpp == 16)
                                        Load_Image_Bitmap16(obj->texture, &bitmap16bit, 0, 0, BITMAP_EXTRACT_MODE_ABS);
                                    else
                                    {
                                        Load_Image_Bitmap(obj->texture, &bitmap16bit, 0, 0, BITMAP_EXTRACT_MODE_ABS);
                                    } // end else 8 bit

                                    // done, so unload the bitmap
                                    Unload_Bitmap_File(&bitmap16bit);

                                    // flag object as having textures
                                    SET_BIT(obj->attr, OBJECT4DV2_ATTR_TEXTURES);

                                } // end if

                                break;
                            } // end if

                        } // end while

                        break;
                    } // end if

                } // end while

                // alright, finally! Now we need to know what the actual shader type, now in the COB format
                // I have decided that in the "reflectance" class that's where we will look at what kind
                // of shader is supposed to be used on the polygon

                //  look for the "Shader class: reflectance"
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader reflectance class not found in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // look for "Shader class: reflectance"
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['class:'] ['reflectance']"))
                    {
                        // now we know the next "shader name" is what we are looking for so, break

                        break;
                    } // end if

                } // end while

                // looking for "Shader name: "xxxxxx" (xxxxxx) " again,
                // and based on it's value we map it to our shading system as follows:
                // "constant" -> MATV1_ATTR_SHADE_MODE_CONSTANT
                // "matte"    -> MATV1_ATTR_SHADE_MODE_FLAT
                // "plastic"  -> MATV1_ATTR_SHADE_MODE_GOURAUD
                // "phong"    -> MATV1_ATTR_SHADE_MODE_FASTPHONG
                // and in the case that in the "color" class, we found a "texture map" then the "shading mode" is
                // "texture map" -> MATV1_ATTR_SHADE_MODE_TEXTURE
                // which must be logically or'ed with the other previous modes
                while (1)
                {
                    // get the next line
                    if (!parser.Getline(PARSER_STRIP_EMPTY_LINES | PARSER_STRIP_WS_ENDS))
                    {
                        //Write_Error("\nshader name ended abruptly! in .COB file %s.", filename);
                        return (0);
                    } // end if

                    // get rid of those quotes
                    ReplaceChars(parser.buffer, parser.buffer, "\"", ' ', 1);

                    // did we find the name?
                    if (parser.Pattern_Match(parser.buffer, "['Shader'] ['name:'] [s>0]"))
                    {
                        // figure out which shader to use
                        if (strcmp(parser.pstrings[2], "constant") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_CONSTANT);
                        } // end if
                        else if (strcmp(parser.pstrings[2], "matte") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_FLAT);
                        } // end if
                        else if (strcmp(parser.pstrings[2], "plastic") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[curr_material + num_materials].attr, MATV1_ATTR_SHADE_MODE_GOURAUD);
                        } // end if
                        else if (strcmp(parser.pstrings[2], "phong") == 0)
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_FASTPHONG);
                        } // end if
                        else
                        {
                            // set the shading mode flag in material
                            SET_BIT(materials[material_index + num_materials].attr, MATV1_ATTR_SHADE_MODE_FLAT);
                        } // end else

                        break;
                    } // end if

                } // end while

                // found the material, break out of while for another pass
                break;

            } // end if found material

        } // end while looking for mat#1

    } // end for curr_material

    // at this point poly_material[] holds all the indices for the polygon materials (zero based, so they need fix up)
    // and we must access the materials array to fill in each polygon with the polygon color, etc.
    // now that we finally have the material libary loaded
    for (int curr_poly = 0; curr_poly < obj->num_polys; curr_poly++)
    {

        // fix up offset
        poly_material[curr_poly] = poly_material[curr_poly] + num_materials;

        // we need to know what color depth we are dealing with, so check
        // the bits per pixel, this assumes that the system has already
        // made the call to DDraw_Init() or set the bit depth
        // if (screen_bpp == 16)
        if (true)
        {
            // cool, 16 bit mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV1_ATTR_RGB16);

            // test if this is a textured poly, if so override the color to WHITE,
            // so we get maximum reflection in lighting stage
            if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_TEXTURE)
                obj->plist[curr_poly].color = RGB16Bit565(255, 255, 255);
            else
                obj->plist[curr_poly].color = RGB16Bit565(materials[poly_material[curr_poly]].color.r,
                                                          materials[poly_material[curr_poly]].color.g,
                                                          materials[poly_material[curr_poly]].color.b);
            //Write_Error("\nPolygon 16-bit");
        } // end
        else
        {
            // 8 bit mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV1_ATTR_8BITCOLOR);

            // test if this is a textured poly, if so override the color to WHITE,
            // so we get maximum reflection in lighting stage
            if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_TEXTURE)
                obj->plist[curr_poly].color = RGBto8BitIndex(255, 255, 255, palette, 0);
            else
                obj->plist[curr_poly].color = RGBto8BitIndex(materials[poly_material[curr_poly]].color.r,
                                                             materials[poly_material[curr_poly]].color.g,
                                                             materials[poly_material[curr_poly]].color.b,
                                                             palette, 0);

            //Write_Error("\nPolygon 8-bit, index=%d", obj->plist[curr_poly].color);
        } // end else

        // now set all the shading flags
        // figure out which shader to use
        if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_CONSTANT)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_CONSTANT);
        } // end if
        else if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_FLAT)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_FLAT);
        } // end if
        else if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_GOURAUD)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_GOURAUD);

            // going to need vertex normals!
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[0]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[1]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[2]].attr, VERTEX4DTV1_ATTR_NORMAL);

        } // end if
        else if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_FASTPHONG)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_FASTPHONG);

            // going to need vertex normals!
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[0]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[1]].attr, VERTEX4DTV1_ATTR_NORMAL);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[2]].attr, VERTEX4DTV1_ATTR_NORMAL);
        } // end if
        else
        {
            // set shading mode to default flat
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_FLAT);

        } // end if

        if (materials[poly_material[curr_poly]].attr & MATV1_ATTR_SHADE_MODE_TEXTURE)
        {
            // set shading mode
            SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_SHADE_MODE_TEXTURE);

            // apply texture to this polygon
            obj->plist[curr_poly].texture = obj->texture;

            // set texture coordinate attributes
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[0]].attr, VERTEX4DTV1_ATTR_TEXTURE);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[1]].attr, VERTEX4DTV1_ATTR_TEXTURE);
            SET_BIT(obj->vlist_local[obj->plist[curr_poly].vert[2]].attr, VERTEX4DTV1_ATTR_TEXTURE);

        } // end if

        // set the material mode to ver. 1.0 emulation (for now only!!!)
        SET_BIT(obj->plist[curr_poly].attr, POLY4DV2_ATTR_DISABLE_MATERIAL);

    } // end for curr_poly

    // local object materials have been added to database, update total materials in system
    num_materials += num_materials_object;

    // now fix up all texture coordinates
    if (obj->texture)
    {
        for (int tvertex = 0; tvertex < num_texture_vertices; tvertex++)
        {
            // step 1: scale the texture coordinates by the texture size
            int texture_size = obj->texture->width;

            // scale 0..1 to 0..texture_size-1
            obj->tlist[tvertex].x *= (texture_size - 1);
            obj->tlist[tvertex].y *= (texture_size - 1);

            // now test for vertex transformation flags
            if (vertex_flags & VERTEX_FLAGS_INVERT_TEXTURE_U)
            {
                obj->tlist[tvertex].x = (texture_size - 1) - obj->tlist[tvertex].x;
            } // end if

            if (vertex_flags & VERTEX_FLAGS_INVERT_TEXTURE_V)
            {
                obj->tlist[tvertex].y = (texture_size - 1) - obj->tlist[tvertex].y;
            } // end if

            if (vertex_flags & VERTEX_FLAGS_INVERT_SWAP_UV)
            {
                float temp;
                SWAP(obj->tlist[tvertex].x, obj->tlist[tvertex].y, temp);
            } // end if

        } // end for

    } // end if there was a texture loaded for this object

#ifdef DEBUG_ON
    for (curr_material = 0; curr_material < num_materials; curr_material++)
    {
        //Write_Error("\nMaterial %d", curr_material);

        //Write_Error("\nint  state    = %d", materials[curr_material].state);
        //Write_Error("\nint  id       = %d", materials[curr_material].id);
        //Write_Error("\nchar name[64] = %s", materials[curr_material].name);
        //Write_Error("\nint  attr     = %d", materials[curr_material].attr);
        //Write_Error("\nint r         = %d", materials[curr_material].color.r);
        //Write_Error("\nint g         = %d", materials[curr_material].color.g);
        //Write_Error("\nint b         = %d", materials[curr_material].color.b);
        //Write_Error("\nint alpha     = %d", materials[curr_material].color.a);
        //Write_Error("\nint color     = %d", materials[curr_material].attr);
        //Write_Error("\nfloat ka      = %f", materials[curr_material].ka);
        //Write_Error("\nkd            = %f", materials[curr_material].kd);
        //Write_Error("\nks            = %f", materials[curr_material].ks);
        //Write_Error("\npower         = %f", materials[curr_material].power);
        //Write_Error("\nchar texture_file = %s\n", materials[curr_material].texture_file);
    } // end for curr_material
#endif

    // now that we know the correct number of polygons, we must allocate memory for them
    // and fix up the object, this is a hack, but the file formats are so stupid by not
    // all starting with NUM_VERTICES, NUM_POLYGONS -- that would make everyone's life
    // easier!

#if 0

// step 1: allocate memory for the polygons
POLY4DV2_PTR plist_temp = NULL;

// allocate memory for polygon list, the correct number of polys was overwritten
// into the object during parsing, so we can use the num_polys field
if (!(plist_temp = (POLY4DV2_PTR)malloc(sizeof(POLY4DV2)*obj->num_polys)))
   return(0);

// step 2:  now copy the polygons into the correct list
memcpy((void *)plist_temp, (void *)obj->plist, sizeof(POLY4DV2));

// step 3: now free the old memory and fix the pointer
free(obj->plist);

// now fix the pointer
obj->plist = plist_temp;

#endif

    // compute the polygon normal lengths
    Compute_OBJECT4DV2_Poly_Normals(obj);

    // compute vertex normals for any gouraud shaded polys
    Compute_OBJECT4DV2_Vertex_Normals(obj);

    // return success
    return (1);

} // end Load_OBJECT4DV2_COB

// CLASS CPARSERV1 IMPLEMENTATION /////////////////////////////////////////

// constructor /////////////////////////////////////////////////

CPARSERV1::CPARSERV1()
{
#ifdef PARSER_DEBUG_ON
    printf("\nEntering CPARSERV1() constructor.");
#endif

    // reset file system
    fstream = NULL;
    Reset();

} // end constructor

// destructor ///////////////////////////////////////////////////

CPARSERV1::~CPARSERV1()
{
#ifdef PARSER_DEBUG_ON
    printf("\nEntering ~CPARSERV1() destructor.");
#endif
    Reset();

} // end destructor

// reset file system ////////////////////////////////////////////
int CPARSERV1::Reset()
{
#ifdef PARSER_DEBUG_ON
    printf("\nEntering Reset().");
#endif

    // reset file buffer
    if (fstream)
        fclose(fstream);

    fstream = NULL;

    // clear and reset buffer
    memset(buffer, 0, sizeof(buffer));
    length = 0;
    num_lines = 0;

    // set comment
    strcpy(comment, PARSER_DEFAULT_COMMENT);

    return (1);

} // end Reset

// open file /////////////////////////////////////////////////////

int CPARSERV1::Open(char *filename)
{
#ifdef PARSER_DEBUG_ON
    printf("\nEntering Open().");
#endif

    // reset file system
    Reset();

    // opens a file
    if ((fstream = fopen(filename, "r")) != NULL)
    {
#ifdef PARSER_DEBUG_ON
        printf("\nOpening file: %s", filename);
#endif
        return (1);
    } // end if
    else
    {
#ifdef PARSER_DEBUG_ON
        printf("\nCouldn't open file: %s", filename);
#endif
        return (0);
    } // end else

} // end Open

// close file ////////////////////////////////////////////////////
int CPARSERV1::Close()
{
    return (Reset());
} // end Close

// get line //////////////////////////////////////////////////////

char *CPARSERV1::Getline(int mode)
{
#ifdef PARSER_DEBUG_ON
    printf("\nEntering Getline().");
#endif

    char *string;

    // gets a single line from the stream
    if (fstream)
    {
        // check translation mode
        if (mode & PARSER_STRIP_EMPTY_LINES)
        {
            // get lines until we get a real one with data on it
            while (1)
            {
                // have we went to the end of the file without getting anything?
                if ((string = fgets(buffer, PARSER_BUFFER_SIZE, fstream)) == NULL)
                    break;

                // we have something, strip ws from it
                int slength = strlen(string);
                int sindex = 0;

                // eat up space
                while (isspace(string[sindex]))
                    sindex++;

                // is there anything left?
                if ((slength - sindex) > 0)
                {
                    // copy the string into place
                    memmove((void *)buffer, (void *)&string[sindex], (slength - sindex) + 1);
                    string = buffer;
                    slength = strlen(string);

                    // strip comments also?
                    if (mode & PARSER_STRIP_COMMENTS)
                    {
                        // does this begin with a comment or end with a comment?
                        char *comment_string = strstr(string, comment);

                        // 3 cases, no comment, comment at beginning, comment at end
                        if (comment_string == NULL)
                            break; // line is valid exit with line

                        // compute index into string from beginning where comment begins
                        int cindex = (int)(comment_string - string);

                        // comment at beginning then continue

                        if (cindex == 0)
                            continue; // this line is a comment, ignore completely, get another
                        else
                        {
                            // comment at end, strip it, insert null where it begins
                            comment_string[0] = 0;
                            break;
                        } // end else

                    } // end if

                    // exit loop, we have something :)
                    break;
                } // end if

            } // end while

        } // end if strip mode
        else
        {
            // just get the next line, don't worry about stripping anything
            string = fgets(buffer, PARSER_BUFFER_SIZE, fstream);
        } // end else

        // was the line valid?
        if (string)
        {
            // increment line count
            num_lines++;

            // final stripping of whitspace
            if (mode & PARSER_STRIP_WS_ENDS)
            {
                StringLtrim(buffer);
                StringRtrim(buffer);
            } // end if

            // compute line length
            length = strlen(buffer);

#ifdef PARSER_DEBUG_ON
            printf("\nString[%d]:%s", length, string);
#endif
            // return the pointer, copy of data already in buffer
            return (string);

        } // end if
        else
        {
#ifdef PARSER_DEBUG_ON
            printf("\nEOF");
#endif
            return (NULL);
        }

    } // end if
    else
    {
#ifdef PARSER_DEBUG_ON
        printf("\nFstream NULL.");
#endif
        return (NULL);
    } // end else

} // end Getline

// sets the comment string ///////////////////////////////////////

int CPARSERV1::SetComment(char *string)
{
    // sets the comment string
    if (strlen(string) < PARSER_MAX_COMMENT)
    {
        strcpy(comment, string);
        return (1);
    } // end if
    else
        return (0);

} // end SetComment

// find pattern in line //////////////////////////////////////////

int CPARSERV1::Pattern_Match(char *string, char *pattern, ...)
{
    // this function tries to match the pattern sent in pattern with
    // the sent string, the results are sent back to the sender in the
    // variable arguments as well as stored in the parameter passing area

    // string literal                        = ['string']
    // floating point number                 = [f]
    // integer number                        = [i]
    // match a string exactly ddd chars      = [s=ddd]
    // match a string less than ddd chars    = [s<ddd]
    // match a string greater than ddd chars = [s>ddd]
    // for example to match "vertex: 34.234 56.34 12.4
    // ['vertex'] [f] [f] [f]

    char token_type[PATTERN_MAX_ARGS];                        // type of token, f,i,s,l
    char token_string[PATTERN_MAX_ARGS][PATTERN_BUFFER_SIZE]; // for literal strings this holds them
    char token_operator[PATTERN_MAX_ARGS];                    // holds and operators for the token, >, <, =, etc.
    int token_numeric[PATTERN_MAX_ARGS];                      // holds any numeric data to qualify the token

    char buffer[PARSER_BUFFER_SIZE]; // working buffer

    // a little error testing
    if ((!string || strlen(string) == 0) || (!pattern || strlen(pattern) == 0))
        return (0);

    // copy line into working area
    strcpy(buffer, string);

    int tok_start = 0,
        tok_end = 0,
        tok_restart = 0,
        tok_first_pass = 0,
        num_tokens = 0;

    // step 1: extract token list
    while (1)
    {
        // eat whitepace
        while (isspace(pattern[tok_start]))
            tok_start++;

        // end of line?
        if (tok_start >= strlen(pattern))
            break;

        // look for beginning of token '['
        if (pattern[tok_start] == '[')
        {
            // now look for token code
            switch (pattern[tok_start + 1])
            {
            case PATTERN_TOKEN_FLOAT: // float
            {
                // make sure token is well formed
                if (pattern[tok_start + 2] != ']')
                    return (0); // error

                // advance token scanner
                tok_start += 3;

                // insert a float into pattern
                token_type[num_tokens] = PATTERN_TOKEN_FLOAT; // type of token, f,i,s,l
                strcpy(token_string[num_tokens], "");         // for literal strings this holds them
                token_operator[num_tokens] = 0;               // holds and operators for the token, >, <, =, etc.
                token_numeric[num_tokens] = 0;

                // increment number of tokens
                num_tokens++;

#ifdef PARSER_DEBUG_ON
                printf("\nFound Float token");
#endif
            }
            break;

            case PATTERN_TOKEN_INT: // integer
            {
                // make sure token is well formed
                if (pattern[tok_start + 2] != ']')
                    return (0); // error

                // advance token scanner
                tok_start += 3;

                // insert a int into pattern
                token_type[num_tokens] = PATTERN_TOKEN_INT; // type of token, f,i,s,l
                strcpy(token_string[num_tokens], "");       // for literal strings this holds them
                token_operator[num_tokens] = 0;             // holds and operators for the token, >, <, =, etc.
                token_numeric[num_tokens] = 0;

                // increment number of tokens
                num_tokens++;

#ifdef PARSER_DEBUG_ON
                printf("\nFound Int token");
#endif
            }
            break;

            case PATTERN_TOKEN_LITERAL: // literal string
            {
                // advance token scanner to begining literal string
                tok_start += 2;
                tok_end = tok_start;

                // eat up string
                while (pattern[tok_end] != PATTERN_TOKEN_LITERAL)
                    tok_end++;

                // make sure string is well formed
                if (pattern[tok_end + 1] != ']')
                    return (0);

                // insert a string into pattern

                // literal string lies from (tok_start - (tok_end-1)
                memcpy(token_string[num_tokens], &pattern[tok_start], (tok_end - tok_start));
                token_string[num_tokens][(tok_end - tok_start)] = 0; // null terminate

                token_type[num_tokens] = PATTERN_TOKEN_LITERAL; // type of token, f,i,s,'
                token_operator[num_tokens] = 0;                 // holds and operators for the token, >, <, =, etc.
                token_numeric[num_tokens] = 0;

#ifdef PARSER_DEBUG_ON
                printf("\nFound Literal token = %s", token_string[num_tokens]);
#endif

                // advance token scanner
                tok_start = tok_end + 2;

                // increment number of tokens
                num_tokens++;
            }
            break;

            case PATTERN_TOKEN_STRING: // ascii string varying length
            {
                // look for comparator
                if (pattern[tok_start + 2] == '=' ||
                    pattern[tok_start + 2] == '>' ||
                    pattern[tok_start + 2] == '<')
                {
                    // extract the number
                    tok_end = tok_start + 3;

                    while (isdigit(pattern[tok_end]))
                        tok_end++;

                    // check for well formed
                    if (pattern[tok_end] != ']')
                        return (0);

                    // copy number in ascii to string and convert to real number
                    memcpy(buffer, &pattern[tok_start + 3], (tok_end - tok_start));
                    buffer[tok_end - tok_start] = 0;

                    // insert a string into pattern
                    token_type[num_tokens] = PATTERN_TOKEN_STRING;       // type of token, f,i,s,l
                    strcpy(token_string[num_tokens], "");                // for literal strings this holds them
                    token_operator[num_tokens] = pattern[tok_start + 2]; // holds and operators for the token, >, <, =, etc.
                    token_numeric[num_tokens] = atoi(buffer);

                } // end if
                else
                    return (0); // not well formed

#ifdef PARSER_DEBUG_ON
                printf("\nFound String token, comparator: %c, characters: %d", token_operator[num_tokens], token_numeric[num_tokens]);
#endif
                // advance token scanner
                tok_start = tok_end + 1;

                // increment number of tokens
                num_tokens++;
            }
            break;

            default:
                break;

            } // end switch

        } // end if

        // end of line?
        if (tok_start >= strlen(pattern))
            break;

    } // end while

#ifdef PARSER_DEBUG_ON
    printf("\nstring to parse: %s", string);
    printf("\nPattern to scan for: %s", pattern);
    printf("\nnumber of tokens found %d", num_tokens);
#endif

    // at this point we have the pattern we need to look for, so look for it
    int pattern_state = PATTERN_STATE_INIT; // initial state for pattern recognizer
    int curr_tok = 0;                       // test for num_tokens
    char token[PATTERN_BUFFER_SIZE];        // token under consideration

    // enter scan state machine
    while (1)
    {
        switch (pattern_state)
        {
        case PATTERN_STATE_INIT:
        {
            // initial state for pattern
            strcpy(buffer, string);

            tok_start = 0;
            tok_end = 0;
            tok_restart = 0;
            tok_first_pass = 1;
            curr_tok = 0;

            // reset output arrays
            num_pints = num_pfloats = num_pstrings = 0;

            // transition to restart
            pattern_state = PATTERN_STATE_RESTART;
        }
        break;

        case PATTERN_STATE_RESTART:
        {
            // pattern may still be here?
            curr_tok = 0;
            tok_first_pass = 1;

            // error detection
            if (tok_end >= strlen(buffer))
                return (0);

            // restart scanner after first token from last pass
            tok_start = tok_end = tok_restart;

            // start validating tokens
            pattern_state = PATTERN_STATE_NEXT;
        }
        break;

        case PATTERN_STATE_NEXT:
        {
            // have we matched pattern yet?
            if (curr_tok >= num_tokens)
            {
                pattern_state = PATTERN_STATE_MATCH;
            }
            else
            {
                // get next token
                if (tok_end >= strlen(buffer))
                    return (0);

                tok_start = tok_end;
                while (isspace(buffer[tok_start]))
                    tok_start++;
                tok_end = tok_start;

                while (!isspace(buffer[tok_end]) && tok_end < strlen(buffer))
                    tok_end++;

                // copy token
                memcpy(token, &buffer[tok_start], tok_end - tok_start);
                token[tok_end - tok_start] = 0;

                // check for error
                if (strlen(token) == 0)
                    return (0);

                // remember position of first token, so we can restart after it on next pass
                // if need
                if (tok_first_pass)
                {
                    tok_first_pass = 0;
                    tok_restart = tok_end;
                } // end if

                // we have the token, set state to check for that token
                switch (token_type[curr_tok])
                {
                case PATTERN_TOKEN_FLOAT:
                {
                    pattern_state = PATTERN_STATE_FLOAT;
                }
                break;
                case PATTERN_TOKEN_INT:
                {
                    pattern_state = PATTERN_STATE_INT;
                }
                break;
                case PATTERN_TOKEN_STRING:
                {
                    pattern_state = PATTERN_STATE_STRING;
                }
                break;
                case PATTERN_TOKEN_LITERAL:
                {
                    pattern_state = PATTERN_STATE_LITERAL;
                }
                break;

                default:
                    break;

                } // end switch

            } // end else
        }
        break;

        case PATTERN_STATE_FLOAT:
        {
            // simply validate this token as a float
            float f = IsFloat(token);

            if (f != FLT_MIN)
            {
                pfloats[num_pfloats++] = f;

                // get next token
                curr_tok++;
                pattern_state = PATTERN_STATE_NEXT;
            } // end if
            else
            {
                // error pattern doesn't match, restart
                pattern_state = PATTERN_STATE_RESTART;
            } // end else
        }
        break;

        case PATTERN_STATE_INT:
        {
            // simply validate this token as a int
            int i = IsInt(token);

            if (i != INT_MIN)
            {
                pints[num_pints++] = i;

                // get next token
                curr_tok++;
                pattern_state = PATTERN_STATE_NEXT;
            } // end if
            else
            {
                // error pattern doesn't match, restart
                pattern_state = PATTERN_STATE_RESTART;
            } // end else
        }
        break;

        case PATTERN_STATE_LITERAL:
        {
            // simply validate this token by comparing to data in table
            if (strcmp(token, token_string[curr_tok]) == 0)
            {
                // increment number of pstrings found and insert into table
                strcpy(pstrings[num_pstrings++], token);

                // get next token
                curr_tok++;
                pattern_state = PATTERN_STATE_NEXT;
            } // end if
            else
            {
                // error pattern doesn't match, restart
                pattern_state = PATTERN_STATE_RESTART;
            } // end else
        }
        break;

        case PATTERN_STATE_STRING:
        {
            // need to test for non-space chars
            // get comparator

            switch (token_operator[curr_tok])
            {
            case '=':
            {
                // we need exactly
                if (strlen(token) == token_numeric[curr_tok])
                {
                    // put this string into table
                    strcpy(pstrings[num_pstrings++], token);

                    // get next token
                    curr_tok++;
                    pattern_state = PATTERN_STATE_NEXT;
                } // end if
                else
                {
                    // error pattern doesn't match, restart
                    pattern_state = PATTERN_STATE_RESTART;
                } // end else
            }
            break;

            case '>':
            {
                // we need greater than
                if (strlen(token) > token_numeric[curr_tok])
                {
                    // put this string into table
                    strcpy(pstrings[num_pstrings++], token);

                    // get next token
                    curr_tok++;
                    pattern_state = PATTERN_STATE_NEXT;
                } // end if
                else
                {
                    // error pattern doesn't match, restart
                    pattern_state = PATTERN_STATE_RESTART;
                } // end else
            }
            break;

            case '<':
            {
                // we need less than
                if (strlen(token) < token_numeric[curr_tok])
                {
                    // put this string into table
                    strcpy(pstrings[num_pstrings++], token);

                    // get next token
                    curr_tok++;
                    pattern_state = PATTERN_STATE_NEXT;
                } // end if
                else
                {
                    // error pattern doesn't match, restart
                    pattern_state = PATTERN_STATE_RESTART;
                } // end else
            }
            break;

            default:
                break;

            } // end switch
        }
        break;

        case PATTERN_STATE_MATCH:
        {
            // we have matched the string, output vars into variable arg list

#ifdef PARSER_DEBUG_ON
            printf("\nPattern: %s matched!", pattern);
#endif

            return (1);
        }
        break;

        case PATTERN_STATE_END:
        {
        }
        break;

        default:
            break;

        } // end switch

    } // end while

} // end Pattern_Match

// END IMPLEMENTATION OF CPARSERV1 CLASS ///////////////////////////////////

void Print_Mat_4X4(MATRIX4X4_PTR ma, char *name = "M")
{
    // prints out a 4x4 matrix
    //Write_Error("\n%s=\n", name);
    //for (int r = 0; r < 4; r++)
    //for (int c = 0; c < 4; c++)
    // Write_Error("%f ", ma->M[r][c]);

} // end Print_Mat_4X4gg

int Init_OBJECT4DV2(OBJECT4DV2_PTR obj, // object to allocate
                    int _num_vertices,
                    int _num_polys,
                    int _num_frames,
                    int destroy)
{
    // this function does nothing more than allocate the memory for an OBJECT4DV2
    // based on the sent data, later we may want to create more robust initializers
    // but the problem is that we don't want to tie the initializer to anthing yet
    // in 99% of cases this all will be done by the call to load object
    // we just might need this function if we manually want to build an object???

    // first destroy the object if it exists
    if (destroy)
        Destroy_OBJECT4DV2(obj);

    // allocate memory for vertex lists
    if (!(obj->vlist_local = (VERTEX4DTV1_PTR)malloc(sizeof(VERTEX4DTV1) * _num_vertices * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->vlist_local, 0, sizeof(VERTEX4DTV1) * _num_vertices * _num_frames);

    if (!(obj->vlist_trans = (VERTEX4DTV1_PTR)malloc(sizeof(VERTEX4DTV1) * _num_vertices * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->vlist_trans, 0, sizeof(VERTEX4DTV1) * _num_vertices * _num_frames);

    // number of texture coordinates always 3*number of polys
    if (!(obj->tlist = (POINT2D_PTR)malloc(sizeof(POINT2D) * _num_polys * 3)))
        return (0);

    // clear data
    memset((void *)obj->tlist, 0, sizeof(POINT2D) * _num_polys * 3);

    // allocate memory for radii arrays
    if (!(obj->avg_radius = (float *)malloc(sizeof(float) * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->avg_radius, 0, sizeof(float) * _num_frames);

    if (!(obj->max_radius = (float *)malloc(sizeof(float) * _num_frames)))
        return (0);

    // clear data
    memset((void *)obj->max_radius, 0, sizeof(float) * _num_frames);

    // allocate memory for polygon list
    if (!(obj->plist = (POLY4DV2_PTR)malloc(sizeof(POLY4DV2) * _num_polys)))
        return (0);

    // clear data
    memset((void *)obj->plist, 0, sizeof(POLY4DV2) * _num_polys);

    // alias head pointers
    obj->head_vlist_local = obj->vlist_local;
    obj->head_vlist_trans = obj->vlist_trans;

    // set some internal variables
    obj->num_frames = _num_frames;
    obj->num_polys = _num_polys;
    obj->num_vertices = _num_vertices;
    obj->total_vertices = _num_vertices * _num_frames;

    // return success
    return (1);

} // end Init_OBJECT4DV2

float Compute_OBJECT4DV2_Radius(OBJECT4DV2_PTR obj)
{
    // this function computes the average and maximum radius for
    // sent object and opdates the object data for the "current frame"
    // it's up to the caller to make sure Set_Frame() for this object
    // has been called to set the object up properly

    // reset incase there's any residue
    obj->avg_radius[obj->curr_frame] = 0;
    obj->max_radius[obj->curr_frame] = 0;

    // loop thru and compute radius
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // update the average and maximum radius
        float dist_to_vertex =
            sqrt(obj->vlist_local[vertex].x * obj->vlist_local[vertex].x +
                 obj->vlist_local[vertex].y * obj->vlist_local[vertex].y +
                 obj->vlist_local[vertex].z * obj->vlist_local[vertex].z);

        // accumulate total radius
        obj->avg_radius[obj->curr_frame] += dist_to_vertex;

        // update maximum radius
        if (dist_to_vertex > obj->max_radius[obj->curr_frame])
            obj->max_radius[obj->curr_frame] = dist_to_vertex;

    } // end for vertex

    // finallize average radius computation
    obj->avg_radius[obj->curr_frame] /= obj->num_vertices;

    // return max radius of frame 0
    return (obj->max_radius[0]);

} // end Compute_OBJECT4DV2_Radius

int ReplaceChars(char *string_in, char *string_out, char *replace_chars, char rep_char, int case_on)
{
    // this function simply replaces the characters from the input string that
    // are listed in replace with the replace char, the results are stored in
    // string_out, string_in and isn't touched, the number of replacments is
    // returned. if case_on = 1 then case is checked, other it's case insensitive

    int num_replacements = 0,            // tracks number of characters replaced
        index_in = 0,                    // curr index into input
        index_out = 0,                   // curr index into output
        sindex,                          // loop var into strip array
        slength = strlen(replace_chars); // length of strip string

    // do some error checking
    if (!string_in || !string_out || strlen(string_in) == 0)
        return (0);

    // nothing to replace
    if (!replace_chars || strlen(replace_chars) == 0)
    {
        strcpy(string_out, string_in);
        return (0);
    } // end if

    // determine if case is important
    if (case_on == 1)
    {
        // perform char by char copy
        while (string_in[index_in])
        {
            for (sindex = 0; sindex < slength; sindex++)
                if (string_in[index_in] == replace_chars[sindex])
                {
                    // replace it
                    string_out[index_out++] = rep_char;
                    index_in++;
                    num_replacements++;
                    break;
                } // end if

            // was a replacement performed?, no just copy then
            if (sindex >= slength)
                string_out[index_out++] = string_in[index_in++];

        } // end while
    }     // end if case_on
    else
    {
        // perform char by char copy with case insensitivity
        while (string_in[index_in])
        {
            for (sindex = 0; sindex < slength; sindex++)
                if (toupper(string_in[index_in]) == toupper(replace_chars[sindex]))
                {
                    // replace it
                    string_out[index_out++] = rep_char;
                    index_in++;
                    num_replacements++;
                    break;
                } // end if

            // was a strip char found?
            if (sindex >= slength)
                string_out[index_out++] = string_in[index_in++];

        } // end while
    }     // end if case_off

    // terminate output string
    string_out[index_out] = 0;

    // return extracts
    return (num_replacements);
}

char *Extract_Filename_From_Path(char *filepath, char *filename)
{
    // this function extracts the filename from a complete path and file
    // "../folder/.../filname.ext"
    // the function operates by scanning backward and looking for the first
    // occurance of "\" or "/" then copies the filename from there to the end
    // test of filepath is valid
    if (!filepath || strlen(filepath) == 0)
        return (NULL);

    int index_end = strlen(filepath) - 1;

    // find filename
    while ((filepath[index_end] != '\\') &&
           (filepath[index_end] != '/') &&
           (filepath[index_end] > 0))
        index_end--;

    // copy file name out into filename var
    memcpy(filename, &filepath[index_end + 1], strlen(filepath) - index_end);

    // return result
    return (filename);

} // end Extract_Filename_From_Path // end ReplaceChars

int Load_Bitmap_File(BITMAP_FILE_PTR bitmap, char *filename)
{
    // this function opens a bitmap file and loads the data into bitmap

    int file_handle, // the file handle
        index;       // looping index

    UCHAR *temp_buffer = NULL; // used to convert 24 bit images to 16 bit
    OFSTRUCT file_data;        // the file data information

    // open the file if it exists
    if ((file_handle = OpenFile(filename, &file_data, OF_READ)) == -1)
        return (0);

    // now load the bitmap file header
    _lread(file_handle, &bitmap->bitmapfileheader, sizeof(BITMAPFILEHEADER));

    // test if this is a bitmap file
    if (bitmap->bitmapfileheader.bfType != BITMAP_ID)
    {
        // close the file
        _lclose(file_handle);

        // return error
        return (0);
    } // end if

    // now we know this is a bitmap, so read in all the sections

    // first the bitmap infoheader

    // now load the bitmap file header
    _lread(file_handle, &bitmap->bitmapinfoheader, sizeof(BITMAPINFOHEADER));

    // now load the color palette if there is one
    if (bitmap->bitmapinfoheader.biBitCount == 8)
    {
        _lread(file_handle, &bitmap->palette, MAX_COLORS_PALETTE * sizeof(PALETTEENTRY));

        // now set all the flags in the palette correctly and fix the reversed
        // BGR RGBQUAD data format
        for (index = 0; index < MAX_COLORS_PALETTE; index++)
        {
            // reverse the red and green fields
            int temp_color = bitmap->palette[index].peRed;
            bitmap->palette[index].peRed = bitmap->palette[index].peBlue;
            bitmap->palette[index].peBlue = temp_color;

            // always set the flags word to this
            bitmap->palette[index].peFlags = PC_NOCOLLAPSE;
        } // end for index

    } // end if

    // finally the image data itself
    _lseek(file_handle, -(int)(bitmap->bitmapinfoheader.biSizeImage), SEEK_END);

    // now read in the image
    if (bitmap->bitmapinfoheader.biBitCount == 8 || bitmap->bitmapinfoheader.biBitCount == 16)
    {
        // delete the last image if there was one
        if (bitmap->buffer)
            free(bitmap->buffer);

        // allocate the memory for the image
        if (!(bitmap->buffer = (UCHAR *)malloc(bitmap->bitmapinfoheader.biSizeImage)))
        {
            // close the file
            _lclose(file_handle);

            // return error
            return (0);
        } // end if

        // now read it in
        _lread(file_handle, bitmap->buffer, bitmap->bitmapinfoheader.biSizeImage);

    } // end if
    else if (bitmap->bitmapinfoheader.biBitCount == 24)
    {
        // allocate temporary buffer to load 24 bit image
        if (!(temp_buffer = (UCHAR *)malloc(bitmap->bitmapinfoheader.biSizeImage)))
        {
            // close the file
            _lclose(file_handle);

            // return error
            return (0);
        } // end if

        // allocate final 16 bit storage buffer
        if (!(bitmap->buffer = (UCHAR *)malloc(2 * bitmap->bitmapinfoheader.biWidth * bitmap->bitmapinfoheader.biHeight)))
        {
            // close the file
            _lclose(file_handle);

            // release working buffer
            free(temp_buffer);

            // return error
            return (0);
        } // end if

        // now read the file in
        _lread(file_handle, temp_buffer, bitmap->bitmapinfoheader.biSizeImage);

        // now convert each 24 bit RGB value into a 16 bit value
        for (index = 0; index < bitmap->bitmapinfoheader.biWidth * bitmap->bitmapinfoheader.biHeight; index++)
        {
            // build up 16 bit color word
            USHORT color;

            // build pixel based on format of directdraw surface
            if (false)
            //    if (dd_pixel_format==DD_PIXEL_FORMAT555)
            {
                // extract RGB components (in BGR order), note the scaling
                UCHAR blue = (temp_buffer[index * 3 + 0] >> 3),
                      green = (temp_buffer[index * 3 + 1] >> 3),
                      red = (temp_buffer[index * 3 + 2] >> 3);
                // use the 555 macro
                //    color = _RGB16BIT555(red,green,blue);
            } // end if 555
            else
                //    if (dd_pixel_format==DD_PIXEL_FORMAT565)
                if (true)
            {
                // extract RGB components (in BGR order), note the scaling
                UCHAR blue = (temp_buffer[index * 3 + 0] >> 3),
                      green = (temp_buffer[index * 3 + 1] >> 2),
                      red = (temp_buffer[index * 3 + 2] >> 3);

                // use the 565 macro
                color = _RGB16BIT565(red, green, blue);

            } // end if 565

            // write color to buffer
            ((USHORT *)bitmap->buffer)[index] = color;

        } // end for index

        // finally write out the correct number of bits
        bitmap->bitmapinfoheader.biBitCount = 16;

        // release working buffer
        free(temp_buffer);

    } // end if 24 bit
    else
    {
        // serious problem
        return (0);

    } // end else

#if 0
// write the file info out 
printf("\nfilename:%s \nsize=%d \nwidth=%d \nheight=%d \nbitsperpixel=%d \ncolors=%d \nimpcolors=%d",
        filename,
        bitmap->bitmapinfoheader.biSizeImage,
        bitmap->bitmapinfoheader.biWidth,
        bitmap->bitmapinfoheader.biHeight,
		bitmap->bitmapinfoheader.biBitCount,
        bitmap->bitmapinfoheader.biClrUsed,
        bitmap->bitmapinfoheader.biClrImportant);
#endif

    // close the file
    _lclose(file_handle);

    // flip the bitmap
    Flip_Bitmap(bitmap->buffer,
                bitmap->bitmapinfoheader.biWidth * (bitmap->bitmapinfoheader.biBitCount / 8),
                bitmap->bitmapinfoheader.biHeight);

    // return success
    return (1);

} // end Load_Bitmap_File

int Create_Bitmap(BITMAP_IMAGE_PTR image, int x, int y, int width, int height, int bpp)
{
    // this function is used to intialize a bitmap, 8 or 16 bit

    // allocate the memory
    if (!(image->buffer = (UCHAR *)malloc(width * height * (bpp >> 3))))
        return (0);

    // initialize variables
    image->state = BITMAP_STATE_ALIVE;
    image->attr = 0;
    image->width = width;
    image->height = height;
    image->bpp = bpp;
    image->x = x;
    image->y = y;
    image->num_bytes = width * height * (bpp >> 3);

    // clear memory out
    memset(image->buffer, 0, width * height * (bpp >> 3));

    // return success
    return (1);

} // end Create_Bitmap

int Load_Image_Bitmap16(BITMAP_IMAGE_PTR image, // bitmap image to load with data
                        BITMAP_FILE_PTR bitmap, // bitmap to scan image data from
                        int cx, int cy,         // cell or absolute pos. to scan image from
                        int mode)               // if 0 then cx,cy is cell position, else
                                                // cx,cy are absolute coords
{
    // this function extracts a 16-bit bitmap out of a 16-bit bitmap file

    // is this a valid bitmap
    if (!image)
        return (0);

    // must be a 16bit bitmap
    USHORT *source_ptr, // working pointers
        *dest_ptr;

    // test the mode of extraction, cell based or absolute
    if (mode == BITMAP_EXTRACT_MODE_CELL)
    {
        // re-compute x,y
        cx = cx * (image->width + 1) + 1;
        cy = cy * (image->height + 1) + 1;
    } // end if

    // extract bitmap data
    source_ptr = (USHORT *)bitmap->buffer +
                 cy * bitmap->bitmapinfoheader.biWidth + cx;

    // assign a pointer to the bimap image
    dest_ptr = (USHORT *)image->buffer;

    int bytes_per_line = image->width * 2;

    // iterate thru each scanline and copy bitmap
    for (int index_y = 0; index_y < image->height; index_y++)
    {
        // copy next line of data to destination
        memcpy(dest_ptr, source_ptr, bytes_per_line);

        // advance pointers
        dest_ptr += image->width;
        source_ptr += bitmap->bitmapinfoheader.biWidth;
    } // end for index_y

    // set state to loaded
    image->attr |= BITMAP_ATTR_LOADED;

    // return success
    return (1);

} // end Load_Image_Bitmap16

int Load_Image_Bitmap(BITMAP_IMAGE_PTR image, // bitmap image to load with data
                      BITMAP_FILE_PTR bitmap, // bitmap to scan image data from
                      int cx, int cy,         // cell or absolute pos. to scan image from
                      int mode)               // if 0 then cx,cy is cell position, else
                                              // cx,cy are absolute coords
{
    // this function extracts a bitmap out of a bitmap file

    // is this a valid bitmap
    if (!image)
        return (0);

    UCHAR *source_ptr, // working pointers
        *dest_ptr;

    // test the mode of extraction, cell based or absolute
    if (mode == BITMAP_EXTRACT_MODE_CELL)
    {
        // re-compute x,y
        cx = cx * (image->width + 1) + 1;
        cy = cy * (image->height + 1) + 1;
    } // end if

    // extract bitmap data
    source_ptr = bitmap->buffer +
                 cy * bitmap->bitmapinfoheader.biWidth + cx;

    // assign a pointer to the bimap image
    dest_ptr = (UCHAR *)image->buffer;

    // iterate thru each scanline and copy bitmap
    for (int index_y = 0; index_y < image->height; index_y++)
    {
        // copy next line of data to destination
        memcpy(dest_ptr, source_ptr, image->width);

        // advance pointers
        dest_ptr += image->width;
        source_ptr += bitmap->bitmapinfoheader.biWidth;
    } // end for index_y

    // set state to loaded
    image->attr |= BITMAP_ATTR_LOADED;

    // return success
    return (1);

} // end Load_Image_Bitmap

int Unload_Bitmap_File(BITMAP_FILE_PTR bitmap)
{
    // this function releases all memory associated with "bitmap"
    if (bitmap->buffer)
    {
        // release memory
        free(bitmap->buffer);

        // reset pointer
        bitmap->buffer = NULL;

    } // end if

    // return success
    return (1);

} // end Unload_Bitmap_File

int RGBto8BitIndex(UCHAR r, UCHAR g, UCHAR b, LPPALETTEENTRY palette, int flush_cache = 0)
{
    // this function hunts thru the loaded 8-bit palette and tries to find the
    // best match to the sent rgb color, the 8-bit index is returned. The algorithm
    // performings a least squares match on the values in the CLUT, also to speed up
    // the process, the last few translated colored are stored in a stack in the format
    // rgbi, so when a new rgb comes in, it is compared against the rgb entries in the
    // table, if found then that index is used, else, the rgb is translated and added to
    // the table, and the table is shifted one slot, so the last element is thrown away,
    // hence the table is FIFO in as much as the first discarded value will be the first
    // this way the system keeps previously translated colors cached, so the fairly long
    // least squared scan doesn't take forever!
    // also note the compression of the RGBI data, a compare is performed on the upper 24 bits only
    // also, if flush_cache = 1 then the local cache is flushed, for example a new palette is loaded

#define COLOR_CACHE_SIZE 16 // 16 entries should do for now

    typedef struct
    {
        UCHAR r, g, b; // the rgb value of this translated color
        UCHAR index;   // the color index that matched is most closely
    } RGBINDEX, *RGBINDEX_PTR;

    static RGBINDEX color_cache[COLOR_CACHE_SIZE]; // the color cache
    static int cache_entries = 0;                  // number of entries in the cache

    // test for flush cache command, new palette coming in...
    if (flush_cache == 1)
        cache_entries = 0;

    // test if the color is in the cache
    for (int cache_index = 0; cache_index < cache_entries; cache_index++)
    {
        // is this a match?
        if (r == color_cache[cache_index].r &&
            g == color_cache[cache_index].g &&
            b == color_cache[cache_index].b)
            return (color_cache[cache_index].index);

    } // end for

    // if we get here then we had no luck, so least sqaures scan for best match
    // and make sure to add results to cache

    int curr_index = -1;       // current color index of best match
    long curr_error = INT_MAX; // distance in color space to nearest match or "error"

    for (int color_index = 0; color_index < 256; color_index++)
    {
        // compute distance to color from target
        long delta_red = abs(palette[color_index].peRed - r);
        long delta_green = abs(palette[color_index].peGreen - g);
        long delta_blue = abs(palette[color_index].peBlue - b);

        long error = (delta_red * delta_red) + (delta_green * delta_green) + (delta_blue * delta_blue);

        // is this color a better match?
        if (error < curr_error)
        {
            curr_index = color_index;
            curr_error = error;
        } // end if

    } // end for color_index

    // at this point we have the new color, insert it into cache
    // shift cache over one entry, copy elements [0 - (n-1)] -> [1 - n]
    memmove((void *)&color_cache[1], (void *)&color_cache[0], COLOR_CACHE_SIZE * sizeof(RGBINDEX) - sizeof(RGBINDEX));

    // now insert the new element
    color_cache[0].r = r;
    color_cache[0].b = b;
    color_cache[0].g = g;
    color_cache[0].index = curr_index;

    // increment number of elements in the cache until saturation
    if (++cache_entries > COLOR_CACHE_SIZE)
        cache_entries = COLOR_CACHE_SIZE;

    // return results
    return (curr_index);

} // end RGBto8BitIndex

int Compute_OBJECT4DV2_Poly_Normals(OBJECT4DV2_PTR obj)
{
    // the normal of a polygon is commonly needed in a number
    // of functions, however, to store a normal turns out to
    // be counterproductive in most cases since the transformation
    // to rotate the normal ends up taking as long as computing the
    // normal -- HOWEVER, if the normal must have unit length, then
    // pre-computing the length of the normal, and then in real-time
    // dividing by this save a length computation, so we get the
    // best of both worlds... thus, this function computes the length
    // of a polygon's normal, but care must be taken, so that we compute
    // the length based on the EXACT same two vectors that all other
    // functions will use when computing the normal
    // in most cases the functions of interest are the lighting functions
    // if we can pre-compute the normal length
    // for all these functions then that will save at least:
    // num_polys_per_frame * (time to compute length of vector)

    // the way we have written the engine, in all cases the normals
    // during lighting are computed as u = v0->v1, and v = v1->v2
    // so as long as we follow that convention we are fine.
    // also, since the new OBJECT4DV2 format supports multiple frames
    // we must perform these calculations for EACH frame of the animation
    // since although the poly indices don't change, the vertice positions
    // do and thus, so do the normals!!!

    // is this object valid
    if (!obj)
        return (0);

    // iterate thru the poly list of the object and compute normals
    // each polygon
    for (int poly = 0; poly < obj->num_polys; poly++)
    {

        // extract vertex indices into master list, rember the polygons are
        // NOT self contained, but based on the vertex list stored in the object
        // itself
        int vindex_0 = obj->plist[poly].vert[0];
        int vindex_1 = obj->plist[poly].vert[1];
        int vindex_2 = obj->plist[poly].vert[2];

        // we need to compute the normal of this polygon face, and recall
        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
        VECTOR4D u, v, n;

        // build u, v
        VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_1].v, &u);
        VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_2].v, &v);

        // compute cross product
        VECTOR4D_Cross(&u, &v, &n);

        // compute length of normal accurately and store in poly nlength
        // +- epsilon later to fix over/underflows
        obj->plist[poly].nlength = VECTOR4D_Length(&n);
    } // end for poly

    // return success
    return (1);

} // end Compute_OBJECT4DV2_Poly_Normals

int Compute_OBJECT4DV2_Vertex_Normals(OBJECT4DV2_PTR obj)
{
    // the vertex normals of each polygon are commonly needed in a number
    // functions, most importantly lighting calculations for gouraud shading
    // however, we only need to compute the vertex normals for polygons that are
    // gouraud shader, so for every vertex we must determine the polygons that
    // share the vertex then compute the average normal, to determine if a polygon
    // contributes we look at the shading flags for the polygon

    // is this object valid
    if (!obj)
        return (0);

    // algorithm: we are going to scan the polygon list and for every polygon
    // that needs normals we are going to "accumulate" the surface normal into all
    // vertices that the polygon touches, and increment a counter to track how many
    // polys contribute to vertex, then when the scan is done the counts will be used
    // to average the accumulated values, so instead of an O(n^2) algorithm, we get a O(c*n)

    // this tracks the polygon indices that touch a particular vertex
    // the array is used to count the number of contributors to the vertex
    // so at the end of the process we can divide each "accumulated" normal
    // and average
    int polys_touch_vertex[OBJECT4DV2_MAX_VERTICES];
    memset((void *)polys_touch_vertex, 0, sizeof(int) * OBJECT4DV2_MAX_VERTICES);

    // iterate thru the poly list of the object, compute its normal, then add
    // each vertice that composes it to the "touching" vertex array
    // while accumulating the normal in the vertex normal array

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        //Write_Error("\nprocessing poly %d", poly);

        // test if this polygon needs vertex normals
        if (obj->plist[poly].attr & POLY4DV2_ATTR_SHADE_MODE_GOURAUD)
        {
            // extract vertex indices into master list, rember the polygons are
            // NOT self contained, but based on the vertex list stored in the object
            // itself
            int vindex_0 = obj->plist[poly].vert[0];
            int vindex_1 = obj->plist[poly].vert[1];
            int vindex_2 = obj->plist[poly].vert[2];

            //Write_Error("\nTouches vertices: %d, %d, %d", vindex_0, vindex_1, vindex_2);

            // we need to compute the normal of this polygon face, and recall
            // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv
            VECTOR4D u, v, n;

            // build u, v
            VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_1].v, &u);
            VECTOR4D_Build(&obj->vlist_local[vindex_0].v, &obj->vlist_local[vindex_2].v, &v);

            // compute cross product
            VECTOR4D_Cross(&u, &v, &n);

            // update vertex array to flag this polygon as a contributor
            polys_touch_vertex[vindex_0]++;
            polys_touch_vertex[vindex_1]++;
            polys_touch_vertex[vindex_2]++;

            //Write_Error("\nPoly touch array v[%d] = %d,  v[%d] = %d,  v[%d] = %d", vindex_0, polys_touch_vertex[vindex_0],                                                                           vindex_1, polys_touch_vertex[vindex_1],                                                                           vindex_2, polys_touch_vertex[vindex_2]);

            // now accumulate the normal into the vertex normal itself
            // note, we do NOT normalize at this point since we want the length of the normal
            // to weight on the average, and since the length is in fact the area of the parallelogram
            // constructed by uxv, so we are taking the "influence" of the area into consideration
            VECTOR4D_Add(&obj->vlist_local[vindex_0].n, &n, &obj->vlist_local[vindex_0].n);
            VECTOR4D_Add(&obj->vlist_local[vindex_1].n, &n, &obj->vlist_local[vindex_1].n);
            VECTOR4D_Add(&obj->vlist_local[vindex_2].n, &n, &obj->vlist_local[vindex_2].n);
        } // end for poly

    } // end if needs vertex normals

    // now we are almost done, we have accumulated all the vertex normals, but need to average them
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        // if this vertex has any contributors then it must need averaging, OR we could check
        // the shading hints flags, they should be one to one
        //Write_Error("\nProcessing vertex: %d, attr: %d, contributors: %d", vertex,                                                                        obj->vlist_local[vertex].attr,                                                                        polys_touch_vertex[vertex]);

        // test if this vertex has a normal and needs averaging
        if (polys_touch_vertex[vertex] >= 1)
        {
            obj->vlist_local[vertex].nx /= polys_touch_vertex[vertex];
            obj->vlist_local[vertex].ny /= polys_touch_vertex[vertex];
            obj->vlist_local[vertex].nz /= polys_touch_vertex[vertex];

            // now normalize the normal
            VECTOR4D_Normalize(&obj->vlist_local[vertex].n);

            //Write_Error("\nAvg Vertex normal: [%f, %f, %f]", obj->vlist_local[vertex].nx,                                                        obj->vlist_local[vertex].ny,                                                        obj->vlist_local[vertex].nz);

        } // end if

    } // end for

    // return success
    return (1);

} // end Compute_OBJECT4DV2_Vertex_Normals

char *StringLtrim(char *string)
{
    // trims whitespace from left side, note is destructive
    int sindex = 0;

    int slength = strlen(string);

    if (!string || slength == 0)
        return (string);

    // trim whitespace by advancing pointer
    while (isspace(string[sindex]) && sindex < slength)
        string[sindex++] = 0; // not needed actually

    // copy string to left
    memmove((void *)string, (void *)&string[sindex], (slength - sindex) + 1);

    // now return pointer
    return (string);

} // end StringLtrim

char *StringRtrim(char *string)
{
    // trims whitespace from right side, note is destructive
    int sindex = 0;

    int slength = strlen(string);

    if (!string || slength == 0)
        return (string);

    // index to end of string
    sindex = slength - 1;

    // trim whitespace by overwriting nulls
    while (isspace(string[sindex]) && sindex >= 0)
        string[sindex--] = 0;

    // string doens't need to be moved, so simply return pointer
    return (string);

} // end StringRtrim

float IsFloat(char *fstring)
{
    // validates the sent string as a float and converts it, if it's not valid
    // the function sends back FLT_MIN, the chances of this being the number
    // validated is slim
    // [whitespace] [sign] [digits] [.digits] [ {d | D | e | E }[sign]digits]

    char *string = fstring;

    // must be of the form
    // [whitespace]
    while (isspace(*string))
        string++;

    // [sign]
    if (*string == '+' || *string == '-')
        string++;

    // [digits]
    while (isdigit(*string))
        string++;

    // [.digits]
    if (*string == '.')
    {
        string++;
        while (isdigit(*string))
            string++;
    }

    // [ {d | D | e | E }[sign]digits]
    if (*string == 'e' || *string == 'E' || *string == 'd' || *string == 'D')
    {
        string++;

        // [sign]
        if (*string == '+' || *string == '-')
            string++;

        // [digits]
        while (isdigit(*string))
            string++;
    }

    // the string better be the same size as the other one
    if (strlen(fstring) == (int)(string - fstring))
        return (atof(fstring));
    else
        return (FLT_MIN);

} // end IsFloat

int IsInt(char *istring)
{
    // validates the sent string as a int and converts it, if it's not valid
    // the function sends back INT_MIN, the chances of this being the number
    // validated is slim
    // [whitespace] [sign]digits

    char *string = istring;

    // must be of the form
    // [whitespace]
    while (isspace(*string))
        string++;

    // [sign]
    if (*string == '+' || *string == '-')
        string++;

    // [digits]
    while (isdigit(*string))
        string++;

    // the string better be the same size as the other one
    if (strlen(istring) == (int)(string - istring))
        return (atoi(istring));
    else
        return (INT_MIN);

} // end IsInt

int Destroy_OBJECT4DV2(OBJECT4DV2_PTR obj) // object to destroy
{
    // this function destroys the sent object, basically frees the memory
    // if any that has been allocated

    // local vertex list
    if (obj->head_vlist_local)
        free(obj->head_vlist_local);

    // transformed vertex list
    if (obj->head_vlist_trans)
        free(obj->head_vlist_trans);

    // texture coordinate list
    if (obj->tlist)
        free(obj->tlist);

    // polygon list
    if (obj->plist)
        free(obj->plist);

    // object radii arrays
    if (obj->avg_radius)
        free(obj->avg_radius);

    if (obj->max_radius)
        free(obj->max_radius);

    // now clear out object completely
    memset((void *)obj, 0, sizeof(OBJECT4DV2));

    // return success
    return (1);

} // end Destroy_OBJECT4DV2

int Flip_Bitmap(UCHAR *image, int bytes_per_line, int height)
{
    // this function is used to flip bottom-up .BMP images

    UCHAR *buffer; // used to perform the image processing
    int index;     // looping index

    // allocate the temporary buffer
    if (!(buffer = (UCHAR *)malloc(bytes_per_line * height)))
        return (0);

    // copy image to work area
    memcpy(buffer, image, bytes_per_line * height);

    // flip vertically
    for (index = 0; index < height; index++)
        memcpy(&image[((height - 1) - index) * bytes_per_line],
               &buffer[index * bytes_per_line], bytes_per_line);

    // release the memory
    free(buffer);

    // return success
    return (1);

} // end Flip_Bitmap

float VECTOR4D_Length(VECTOR4D_PTR va)
{
    // computes the magnitude of a vector, slow

    return (sqrtf(va->x * va->x + va->y * va->y + va->z * va->z));

} // end VECTOR4D_Length

void Reset_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list)
{
    rend_list->num_polys = 0; // that was hard!
} // end Reset_RENDERLIST4DV2

void Reset_OBJECT4DV2(OBJECT4DV2_PTR obj)
{
    // this function resets the sent object and redies it for
    // transformations, basically just resets the culled, clipped and
    // backface flags, but here's where you would add stuff
    // to ready any object for the pipeline
    // the object is valid, let's rip it apart polygon by polygon
    // note: works on the entire object, all frames

    // reset object's culled flag
    RESET_BIT(obj->state, OBJECT4DV2_STATE_CULLED);

    // now the clipped and backface flags for the polygons
    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        // acquire polygon
        POLY4DV2_PTR curr_poly = &obj->plist[poly];

        // first is this polygon even visible?
        if (!(curr_poly->state & POLY4DV2_STATE_ACTIVE))
            continue; // move onto next poly

        // reset clipped and backface flags
        RESET_BIT(curr_poly->state, POLY4DV2_STATE_CLIPPED);
        RESET_BIT(curr_poly->state, POLY4DV2_STATE_BACKFACE);
        RESET_BIT(curr_poly->state, POLY4DV2_STATE_LIT);

    } // end for poly

} // end Reset_OBJECT4DV2

void Transform_OBJECT4DV2(OBJECT4DV2_PTR obj,  // object to transform
                          MATRIX4X4_PTR mt,    // transformation matrix
                          int coord_select,    // selects coords to transform
                          int transform_basis, // flags if vector orientation
                                               // should be transformed too
                          int all_frames)      // should all frames be transformed

{
    // this function simply transforms all of the vertices in the local or trans
    // array by the sent matrix, since the object may have multiple frames, it
    // takes that into consideration
    // also vertex normals are rotated, however, if there is a translation factor
    // in the sent matrix that will corrupt the normals, later we might want to
    // null out the last row of the matrix before transforming the normals?
    // future optimization: set flag in object attributes, and objects without
    // vertex normals can be rotated without the test in line

    // single frame or all frames?
    if (!all_frames)
    {
        // what coordinates should be transformed?
        switch (coord_select)
        {
        case TRANSFORM_LOCAL_ONLY:
        {
            // transform each local/model vertex of the object mesh in place
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->vlist_local[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->vlist_local[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_TRANS_ONLY:
        {
            // transform each "transformed" vertex of the object mesh in place
            // remember, the idea of the vlist_trans[] array is to accumulate
            // transformations
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->vlist_trans[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->vlist_trans[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->vlist_trans[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_LOCAL_TO_TRANS:
        {
            // transform each local/model vertex of the object mesh and store result
            // in "transformed" vertex list
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].v, mt, &obj->vlist_trans[vertex].v);

                // transform vertex normal if needed
                if (obj->vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform point
                    Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex].n, mt, &obj->vlist_trans[vertex].n);
                } // end if

            } // end for index
        }
        break;

        default:
            break;

        } // end switch

    }    // end if single frame
    else // transform all frames
    {
        // what coordinates should be transformed?
        switch (coord_select)
        {
        case TRANSFORM_LOCAL_ONLY:
        {
            // transform each local/model vertex of the object mesh in place
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->head_vlist_local[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->head_vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->head_vlist_local[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_TRANS_ONLY:
        {
            // transform each "transformed" vertex of the object mesh in place
            // remember, the idea of the vlist_trans[] array is to accumulate
            // transformations
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_trans[vertex].v, mt, &presult);

                // store result back
                VECTOR4D_COPY(&obj->head_vlist_trans[vertex].v, &presult);

                // transform vertex normal if needed
                if (obj->head_vlist_trans[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform normal
                    Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_trans[vertex].n, mt, &presult);

                    // store result back
                    VECTOR4D_COPY(&obj->head_vlist_trans[vertex].n, &presult);
                } // end if

            } // end for index
        }
        break;

        case TRANSFORM_LOCAL_TO_TRANS:
        {
            // transform each local/model vertex of the object mesh and store result
            // in "transformed" vertex list
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].v, mt, &obj->head_vlist_trans[vertex].v);

                // transform vertex normal if needed
                if (obj->head_vlist_local[vertex].attr & VERTEX4DTV1_ATTR_NORMAL)
                {
                    // transform point
                    Mat_Mul_VECTOR4D_4X4(&obj->head_vlist_local[vertex].n, mt, &obj->head_vlist_trans[vertex].n);
                } // end if

            } // end for index
        }
        break;

        default:
            break;

        } // end switch

    } // end else multiple frames

    // finally, test if transform should be applied to orientation basis
    // hopefully this is a rotation, otherwise the basis will get corrupted
    if (transform_basis)
    {
        // now rotate orientation basis for object
        VECTOR4D vresult; // use to rotate each orientation vector axis

        // rotate ux of basis
        Mat_Mul_VECTOR4D_4X4(&obj->ux, mt, &vresult);
        VECTOR4D_COPY(&obj->ux, &vresult);

        // rotate uy of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uy, mt, &vresult);
        VECTOR4D_COPY(&obj->uy, &vresult);

        // rotate uz of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uz, mt, &vresult);
        VECTOR4D_COPY(&obj->uz, &vresult);
    } // end if

} // end Transform_OBJECT4DV2

void Model_To_World_OBJECT4DV2(OBJECT4DV2_PTR obj,
                               int coord_select,
                               int all_frames)
{
    // NOTE: Not matrix based
    // this function converts the local model coordinates of the
    // sent object into world coordinates, the results are stored
    // in the transformed vertex list (vlist_trans) within the object

    // interate thru vertex list and transform all the model/local
    // coords to world coords by translating the vertex list by
    // the amount world_pos and storing the results in vlist_trans[]
    // no need to transform vertex normals, they are invariant of position

    if (!all_frames)
    {
        if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
        {
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&obj->vlist_local[vertex].v, &obj->world_pos, &obj->vlist_trans[vertex].v);
                // copy normal
                VECTOR4D_COPY(&obj->vlist_trans[vertex].n, &obj->vlist_local[vertex].n);

            } // end for vertex
        }     // end if local
        else
        { // TRANSFORM_TRANS_ONLY
            for (int vertex = 0; vertex < obj->num_vertices; vertex++)
            {
                // std::cout<<"TRANSFORM_TRANS_ONLY"<<std::endl;
                // translate vertex
                VECTOR4D_Add(&obj->vlist_trans[vertex].v, &obj->world_pos, &obj->vlist_trans[vertex].v);
            } // end for vertex
        }     // end else trans

    }    // end if single frame
    else // all frames
    {
        if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
        {
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&obj->head_vlist_local[vertex].v, &obj->world_pos, &obj->head_vlist_trans[vertex].v);
                // copy normal
                VECTOR4D_COPY(&obj->head_vlist_trans[vertex].n, &obj->head_vlist_local[vertex].n);
            } // end for vertex
        }     // end if local
        else
        { // TRANSFORM_TRANS_ONLY
            for (int vertex = 0; vertex < obj->total_vertices; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&obj->head_vlist_trans[vertex].v, &obj->world_pos, &obj->head_vlist_trans[vertex].v);
            } // end for vertex
        }     // end else trans

    } // end if all frames

} // end Model_To_World_OBJECT4DV2

int Insert_OBJECT4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                     OBJECT4DV2_PTR obj,
                                     int insert_local = 0)

{
    // { andre work in progress, rewrite with materials...}

    // converts the entire object into a face list and then inserts
    // the visible, active, non-clipped, non-culled polygons into
    // the render list, also note the flag insert_local control
    // whether or not the vlist_local or vlist_trans vertex list
    // is used, thus you can insert an object "raw" totally untranformed
    // if you set insert_local to 1, default is 0, that is you would
    // only insert an object after at least the local to world transform
    // the last parameter is used to control if their has been
    // a lighting step that has generated a light value stored
    // in the upper 16-bits of color, if lighting_on = 1 then
    // this value is used to overwrite the base color of the
    // polygon when its sent to the rendering list

    unsigned int base_color; // save base color of polygon

    // is this objective inactive or culled or invisible?
    if (!(obj->state & OBJECT4DV2_STATE_ACTIVE) ||
        (obj->state & OBJECT4DV2_STATE_CULLED) ||
        !(obj->state & OBJECT4DV2_STATE_VISIBLE))
        return (0);

    // the object is valid, let's rip it apart polygon by polygon
    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        // acquire polygon
        POLY4DV2_PTR curr_poly = &obj->plist[poly];
        // std::cout << curr_poly->color << std::endl;

        // first is this polygon even visible?
        if (!(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // override vertex list polygon refers to
        // the case that you want the local coords used
        // first save old pointer
        VERTEX4DTV1_PTR vlist_old = curr_poly->vlist;

        if (insert_local)
            curr_poly->vlist = obj->vlist_local;
        else
            curr_poly->vlist = obj->vlist_trans;

        // 丢到了 obj->vlist_trans[vertex].v
        // now insert this polygon
        if (!Insert_POLY4DV2_RENDERLIST4DV2(rend_list, curr_poly))
        {
            // fix vertex list pointer
            curr_poly->vlist = vlist_old;

            // the whole object didn't fit!
            return (0);
        } // end if

        // fix vertex list pointer
        curr_poly->vlist = vlist_old;

    } // end for

    // return success
    return (1);

} // end Insert_OBJECT4DV2_RENDERLIST4DV2

int Insert_POLY4DV2_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                   POLY4DV2_PTR poly)
{
    // converts the sent POLY4DV2 into a POLYF4DV2 and inserts it
    // into the render list, this function needs optmizing

    // step 0: are we full?
    if (rend_list->num_polys >= RENDERLIST4DV2_MAX_POLYS)
        return (0);

    // step 1: copy polygon into next opening in polygon render list

    // point pointer to polygon structure
    rend_list->poly_ptrs[rend_list->num_polys] = &rend_list->poly_data[rend_list->num_polys];

    // copy fields { ??????????? make sure ALL fields are copied, normals, textures, etc!!!  }
    rend_list->poly_data[rend_list->num_polys].state = poly->state;
    rend_list->poly_data[rend_list->num_polys].attr = poly->attr;
    rend_list->poly_data[rend_list->num_polys].color = poly->color;
    rend_list->poly_data[rend_list->num_polys].nlength = poly->nlength;
    rend_list->poly_data[rend_list->num_polys].texture = poly->texture;

    // poly could be lit, so copy these too...
    rend_list->poly_data[rend_list->num_polys].lit_color[0] = poly->lit_color[0];
    rend_list->poly_data[rend_list->num_polys].lit_color[1] = poly->lit_color[1];
    rend_list->poly_data[rend_list->num_polys].lit_color[2] = poly->lit_color[2];

    // now copy vertices, be careful! later put a loop, but for now
    // know there are 3 vertices always!
    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[0],
                     &poly->vlist[poly->vert[0]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[1],
                     &poly->vlist[poly->vert[1]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].tvlist[2],
                     &poly->vlist[poly->vert[2]]);

    // and copy into local vertices too
    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[0],
                     &poly->vlist[poly->vert[0]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[1],
                     &poly->vlist[poly->vert[1]]);

    VERTEX4DTV1_COPY(&rend_list->poly_data[rend_list->num_polys].vlist[2],
                     &poly->vlist[poly->vert[2]]);

    // finally the texture coordinates, this has to be performed manually
    // since at this point in the pipeline the vertices do NOT have texture
    // coordinate, the polygons DO, however, now, there are 3 vertices for
    // EVERY polygon, rather than vertex sharing, so we can copy the texture
    // coordinates out of the indexed arrays into the VERTEX4DTV1 structures
    rend_list->poly_data[rend_list->num_polys].tvlist[0].t = poly->tlist[poly->text[0]];
    rend_list->poly_data[rend_list->num_polys].tvlist[1].t = poly->tlist[poly->text[1]];
    rend_list->poly_data[rend_list->num_polys].tvlist[2].t = poly->tlist[poly->text[2]];

    rend_list->poly_data[rend_list->num_polys].vlist[0].t = poly->tlist[poly->text[0]];
    rend_list->poly_data[rend_list->num_polys].vlist[1].t = poly->tlist[poly->text[1]];
    rend_list->poly_data[rend_list->num_polys].vlist[2].t = poly->tlist[poly->text[2]];

    // now the polygon is loaded into the next free array position, but
    // we need to fix up the links

    // test if this is the first entry
    if (rend_list->num_polys == 0)
    {
        // set pointers to null, could loop them around though to self
        rend_list->poly_data[0].next = NULL;
        rend_list->poly_data[0].prev = NULL;
    } // end if
    else
    {
        // first set this node to point to previous node and next node (null)
        rend_list->poly_data[rend_list->num_polys].next = NULL;
        rend_list->poly_data[rend_list->num_polys].prev =
            &rend_list->poly_data[rend_list->num_polys - 1];

        // now set previous node to point to this node
        rend_list->poly_data[rend_list->num_polys - 1].next =
            &rend_list->poly_data[rend_list->num_polys];
    } // end else

    // increment number of polys in list
    rend_list->num_polys++;

    // return successful insertion
    return (1);

} // end Insert_POLY4DV2_RENDERLIST4DV2

void Remove_Backfaces_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function removes the backfaces from polygon list
    // the function does this based on the polygon list data
    // tvlist along with the camera position (only)
    // note that only the backface state is set in each polygon

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // test this polygon if and only if it's not clipped, not culled,
        // active, and visible and not 2 sided. Note we test for backface in the event that
        // a previous call might have already determined this, so why work
        // harder!
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->attr & POLY4DV2_ATTR_2SIDED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // we need to compute the normal of this polygon face, and recall
        // that the vertices are in cw order, u = p0->p1, v=p0->p2, n=uxv
        VECTOR4D u, v, n;

        // build u, v
        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

        // compute cross product
        VECTOR4D_Cross(&u, &v, &n);

        // now create eye vector to viewpoint
        VECTOR4D view;
        VECTOR4D_Build(&curr_poly->tvlist[0].v, &cam->pos, &view);

        // and finally, compute the dot product
        float dp = VECTOR4D_Dot(&n, &view);

        // if the sign is > 0 then visible, 0 = scathing, < 0 invisible
        if (dp <= 0.0)
            SET_BIT(curr_poly->state, POLY4DV2_STATE_BACKFACE);

    } // end for poly

} // end Remove_Backfaces_RENDERLIST4DV2

//计算出当前poly的color
int Light_RENDERLIST4DV2_World16(RENDERLIST4DV2_PTR rend_list, // list to process
                                 CAM4DV1_PTR cam,              // camera position
                                 LIGHTV1_PTR lights,           // light list (might have more than one)
                                 int max_lights)               // maximum lights in list
{

    // 16-bit version of function
    // function lights the entire rendering list based on the sent lights and camera. the function supports
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
    // this function now performs emissive, flat, and gouraud lighting, results are stored in the
    // lit_color[] array of each polygon

    unsigned int r_base, g_base, b_base, // base color being lit
        r_sum, g_sum, b_sum,             // sum of lighting process over all lights
        r_sum0, g_sum0, b_sum0,
        r_sum1, g_sum1, b_sum1,
        r_sum2, g_sum2, b_sum2,
        ri, gi, bi,
        shaded_color; // final color

    float dp, // dot product
        dist, // distance from light to surface
        dists,
        i,     // general intensities
        nl,    // length of normal
        atten; // attenuation computations

    VECTOR4D u, v, n, l, d, s; // used for cross product and light vector calculations

    //Write_Error("\nEntering lighting function");

    // for each valid poly, light it...
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // light this polygon if and only if it's not clipped, not culled,
        // active, and visible
        if (!(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE) ||
            (curr_poly->state & POLY4DV2_STATE_LIT))
            continue; // move onto next poly

            //Write_Error("\npoly %d",poly);

#ifdef DEBUG_ON
        // track rendering stats
        debug_polys_lit_per_frame++;
#endif

        // set state of polygon to lit
        SET_BIT(curr_poly->state, POLY4DV2_STATE_LIT);

        // we will use the transformed polygon vertex list since the backface removal
        // only makes sense at the world coord stage further of the pipeline

        // test the lighting mode of the polygon (use flat for flat, gouraud))
        if (curr_poly->attr & POLY4DV2_ATTR_SHADE_MODE_FLAT)
        {
            //Write_Error("\nEntering Flat Shader");

            // step 1: extract the base color out in RGB mode
            // assume 565 format
            _RGB565FROM16BIT(curr_poly->color, &r_base, &g_base, &b_base);
            // scale to 8 bit
            r_base <<= 3;
            g_base <<= 2;
            b_base <<= 3;

            // std::cout<<"color = " << r_base<<" "<< g_base <<" "<< b_base<<std::endl;
            //Write_Error("\nBase color=%d,%d,%d", r_base, g_base, b_base);

            // initialize color sum
            r_sum = 0;
            g_sum = 0;
            b_sum = 0;

            //Write_Error("\nsum color=%d,%d,%d", r_sum, g_sum, b_sum);

            // new optimization:
            // when there are multiple lights in the system we will end up performing numerous
            // redundant calculations to minimize this my strategy is to set key variables to
            // to MAX values on each loop, then during the lighting calcs to test the vars for
            // the max value, if they are the max value then the first light that needs the math
            // will do it, and then save the information into the variable (causing it to change state
            // from an invalid number) then any other lights that need the math can use the previously
            // computed value

            // set surface normal.z to FLT_MAX to flag it as non-computed
            n.z = FLT_MAX;

            // loop thru lights
            for (int curr_light = 0; curr_light < max_lights; curr_light++)
            {
                // is this light active
                if (lights[curr_light].state == LIGHTV1_STATE_OFF)
                    continue;

                //Write_Error("\nprocessing light %d",curr_light);

                // what kind of light are we dealing with
                if (lights[curr_light].attr & LIGHTV1_ATTR_AMBIENT)
                {
                    //Write_Error("\nEntering ambient light...");

                    // simply multiply each channel against the color of the
                    // polygon then divide by 256 to scale back to 0..255
                    // use a shift in real life!!! >> 8
                    r_sum += ((lights[curr_light].c_ambient.r * r_base) / 256);
                    g_sum += ((lights[curr_light].c_ambient.g * g_base) / 256);
                    b_sum += ((lights[curr_light].c_ambient.b * b_base) / 256);

                    //Write_Error("\nambient sum=%d,%d,%d", r_sum, g_sum, b_sum);

                    // there better only be one ambient light!

                }                                                         // end if
                else if (lights[curr_light].attr & LIGHTV1_ATTR_INFINITE) ///////////////////////////////////////////
                {
                    //Write_Error("\nEntering infinite light...");

                    // infinite lighting, we need the surface normal, and the direction
                    // of the light source

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

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

                    //Write_Error("\ninfinite sum=%d,%d,%d", r_sum, g_sum, b_sum);

                }                                                      // end if infinite light
                else if (lights[curr_light].attr & LIGHTV1_ATTR_POINT) ///////////////////////////////////////
                {
                    //Write_Error("\nEntering point light...");

                    // perform point light computations
                    // light model for point light is once again:
                    //              I0point * Clpoint
                    //  I(d)point = ___________________
                    //              kc +  kl*d + kq*d2
                    //
                    //  Where d = |p - s|
                    // thus it's almost identical to the infinite light, but attenuates as a function
                    // of distance from the point source to the surface point being lit

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

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

                    //Write_Error("\npoint sum=%d,%d,%d",r_sum,g_sum,b_sum);

                }                                                           // end if point
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT1) ////////////////////////////////////
                {
                    //Write_Error("\nentering spot light1...");

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

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

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

                    //Write_Error("\nspotlight sum=%d,%d,%d",r_sum, g_sum, b_sum);

                }                                                           // end if spotlight1
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT2) // simple version ////////////////////
                {
                    //Write_Error("\nEntering spotlight2 ...");

                    // perform spot light computations
                    // light model for spot light simple version is once again:
                    //         	     I0spotlight * Clspotlight * MAX( (l . s), 0)^pf
                    // I(d)spotlight = __________________________________________
                    //               		 kc + kl*d + kq*d2
                    // Where d = |p - s|, and pf = power factor

                    // thus it's almost identical to the point, but has the extra term in the numerator
                    // relating the angle between the light source and the point on the surface

                    // test if we already computed poly normal in previous calculation
                    if (n.z == FLT_MAX)
                    {
                        // we need to compute the normal of this polygon face, and recall
                        // that the vertices are in cw order, u=p0->p1, v=p0->p2, n=uxv

                        // build u, v
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[1].v, &u);
                        VECTOR4D_Build(&curr_poly->tvlist[0].v, &curr_poly->tvlist[2].v, &v);

                        // compute cross product
                        VECTOR4D_Cross(&u, &v, &n);
                    } // end if

                    // at this point, we are almost ready, but we have to normalize the normal vector!
                    // this is a key optimization we can make later, we can pre-compute the length of all polygon
                    // normals, so this step can be optimized
                    // compute length of normal
                    //nl = VECTOR4D_Length_Fast2(&n);
                    nl = curr_poly->nlength;

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
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[0].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

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

                    //Write_Error("\nSpotlight sum=%d,%d,%d",r_sum, g_sum, b_sum);

                } // end if spot light

            } // end for light

            // make sure colors aren't out of range
            if (r_sum > 255)
                r_sum = 255;
            if (g_sum > 255)
                g_sum = 255;
            if (b_sum > 255)
                b_sum = 255;

            //Write_Error("\nWriting final values to polygon %d = %d,%d,%d", poly, r_sum, g_sum, b_sum);

            // write the color over current color
            curr_poly->lit_color[0] = RGB16Bit565(r_sum, g_sum, b_sum);

        }                                                            // end if
        else if (curr_poly->attr & POLY4DV2_ATTR_SHADE_MODE_GOURAUD) /////////////////////////////////
        {
            // gouraud shade, unfortunetly at this point in the pipeline, we have lost the original
            // mesh, and only have triangles, thus, many triangles will share the same vertices and
            // they will get lit 2x since we don't have any way to tell this, alas, performing lighting
            // at the object level is a better idea when gouraud shading is performed since the
            // commonality of vertices is still intact, in any case, lighting here is similar to polygon
            // flat shaded, but we do it 3 times, once for each vertex, additionally there are lots
            // of opportunities for optimization, but I am going to lay off them for now, so the code
            // is intelligible, later we will optimize

            //Write_Error("\nEntering gouraud shader...");

            // step 1: extract the base color out in RGB mode
            // assume 565 format
            _RGB565FROM16BIT(curr_poly->color, &r_base, &g_base, &b_base);

            // scale to 8 bit
            r_base <<= 3;
            g_base <<= 2;
            b_base <<= 3;

            //Write_Error("\nBase color=%d, %d, %d", r_base, g_base, b_base);

            // initialize color sum(s) for vertices
            r_sum0 = 0;
            g_sum0 = 0;
            b_sum0 = 0;

            r_sum1 = 0;
            g_sum1 = 0;
            b_sum1 = 0;

            r_sum2 = 0;
            g_sum2 = 0;
            b_sum2 = 0;

            //Write_Error("\nColor sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,   r_sum1, g_sum1, b_sum1, r_sum2, g_sum2, b_sum2);

            // new optimization:
            // when there are multiple lights in the system we will end up performing numerous
            // redundant calculations to minimize this my strategy is to set key variables to
            // to MAX values on each loop, then during the lighting calcs to test the vars for
            // the max value, if they are the max value then the first light that needs the math
            // will do it, and then save the information into the variable (causing it to change state
            // from an invalid number) then any other lights that need the math can use the previously
            // computed value, however, since we already have the normals, not much here to cache on
            // a large scale, but small scale stuff is there, however, we will optimize those later

            // loop thru lights
            for (int curr_light = 0; curr_light < max_lights; curr_light++)
            {
                // is this light active
                if (lights[curr_light].state == LIGHTV1_STATE_OFF)
                    continue;

                //Write_Error("\nprocessing light %d", curr_light);

                // what kind of light are we dealing with
                if (lights[curr_light].attr & LIGHTV1_ATTR_AMBIENT) ///////////////////////////////
                {
                    //Write_Error("\nEntering ambient light....");

                    // simply multiply each channel against the color of the
                    // polygon then divide by 256 to scale back to 0..255
                    // use a shift in real life!!! >> 8
                    ri = ((lights[curr_light].c_ambient.r * r_base) / 256);
                    gi = ((lights[curr_light].c_ambient.g * g_base) / 256);
                    bi = ((lights[curr_light].c_ambient.b * b_base) / 256);

                    // ambient light has the same affect on each vertex
                    r_sum0 += ri;
                    g_sum0 += gi;
                    b_sum0 += bi;

                    r_sum1 += ri;
                    g_sum1 += gi;
                    b_sum1 += bi;

                    r_sum2 += ri;
                    g_sum2 += gi;
                    b_sum2 += bi;

                    // there better only be one ambient light!
                    //Write_Error("\nexiting ambient ,sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                }                                                         // end if
                else if (lights[curr_light].attr & LIGHTV1_ATTR_INFINITE) /////////////////////////////////
                {
                    //Write_Error("\nentering infinite light...");

                    // infinite lighting, we need the surface normal, and the direction
                    // of the light source

                    // no longer need to compute normal or length, we already have the vertex normal
                    // and it's length is 1.0
                    // ....

                    // ok, recalling the lighting model for infinite lights
                    // I(d)dir = I0dir * Cldir
                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // need to perform lighting for each vertex (lots of redundant math, optimize later!)

                    //Write_Error("\nv0=[%f, %f, %f]=%f, v1=[%f, %f, %f]=%f, v2=[%f, %f, %f]=%f",
                    // curr_poly->tvlist[0].n.x, curr_poly->tvlist[0].n.y,curr_poly->tvlist[0].n.z, VECTOR4D_Length(&curr_poly->tvlist[0].n),
                    // curr_poly->tvlist[1].n.x, curr_poly->tvlist[1].n.y,curr_poly->tvlist[1].n.z, VECTOR4D_Length(&curr_poly->tvlist[1].n),
                    // curr_poly->tvlist[2].n.x, curr_poly->tvlist[2].n.y,curr_poly->tvlist[2].n.z, VECTOR4D_Length(&curr_poly->tvlist[2].n) );

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp;
                        r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp;
                        r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        i = 128 * dp;
                        r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    //Write_Error("\nexiting infinite, color sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                }                                                      // end if infinite light
                else if (lights[curr_light].attr & LIGHTV1_ATTR_POINT) //////////////////////////////////////
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

                    // .. normal already in vertex

                    //Write_Error("\nEntering point light....");

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // perform the calculation for all 3 vertices

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (dist * atten);


                        r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        // r_sum0 += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                        // g_sum0 += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                        // b_sum0 += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (dist * atten);


                        r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        // r_sum1 += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                        // g_sum1 += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                        // b_sum1 += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
                    } // end if

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &l);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (dist * atten);
                        r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        // r_sum2 += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                        // g_sum2 += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                        // b_sum2 += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
                    } // end if

                    //Write_Error("\nexiting point light, rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1, r_sum2, g_sum2, b_sum2);

                }                                                           // end if point
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT1) ///////////////////////////////////////
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

                    //Write_Error("\nentering spotlight1....");

                    // .. normal is already computed

                    // compute vector from surface to light
                    VECTOR4D_Build(&curr_poly->tvlist[0].v, &lights[curr_light].pos, &l);

                    // compute distance and attenuation
                    dist = VECTOR4D_Length_Fast2(&l);

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    // note that I use the direction of the light here rather than a the vector to the light
                    // thus we are taking orientation into account which is similar to the spotlight model

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (atten);

                        r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (atten);

                        r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end i

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        atten = (lights[curr_light].kc + lights[curr_light].kl * dist + lights[curr_light].kq * dist * dist);

                        i = 128 * dp / (atten);

                        r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                        g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                        b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                    } // end i

                    //Write_Error("\nexiting spotlight1, sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                }                                                           // end if spotlight1
                else if (lights[curr_light].attr & LIGHTV1_ATTR_SPOTLIGHT2) // simple version //////////////////////////
                {
                    // perform spot light computations
                    // light model for spot light simple version is once again:
                    //         	     I0spotlight * Clspotlight * MAX( (l . s), 0)^pf
                    // I(d)spotlight = __________________________________________
                    //               		 kc + kl*d + kq*d2
                    // Where d = |p - s|, and pf = power factor

                    // thus it's almost identical to the point, but has the extra term in the numerator
                    // relating the angle between the light source and the point on the surface

                    // .. already have normals and length are 1.0

                    // and for the diffuse model
                    // Itotald =   Rsdiffuse*Idiffuse * (n . l)
                    // so we basically need to multiple it all together
                    // notice the scaling by 128, I want to avoid floating point calculations, not because they
                    // are slower, but the conversion to and from cost cycles

                    //Write_Error("\nEntering spotlight2...");

                    // tons of redundant math here! lots to optimize later!

                    // vertex 0
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[0].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[0].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (atten);

                            r_sum0 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum0 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum0 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);

                        } // end if

                    } // end if

                    // vertex 1
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[1].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[1].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (atten);

                            r_sum1 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum1 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum1 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);

                        } // end if

                    } // end if

                    // vertex 2
                    dp = VECTOR4D_Dot(&curr_poly->tvlist[2].n, &lights[curr_light].dir);

                    // only add light if dp > 0
                    if (dp > 0)
                    {
                        // compute vector from light to surface (different from l which IS the light dir)
                        VECTOR4D_Build(&lights[curr_light].pos, &curr_poly->tvlist[2].v, &s);

                        // compute length of s (distance to light source) to normalize s for lighting calc
                        dists = VECTOR4D_Length_Fast2(&s);

                        // compute spot light term (s . l)
                        float dpsl = VECTOR4D_Dot(&s, &lights[curr_light].dir) / dists;

                        // proceed only if term is positive
                        if (dpsl > 0)
                        {
                            // compute attenuation
                            atten = (lights[curr_light].kc + lights[curr_light].kl * dists + lights[curr_light].kq * dists * dists);

                            // for speed reasons, pf exponents that are less that 1.0 are out of the question, and exponents
                            // must be integral
                            float dpsl_exp = dpsl;

                            // exponentiate for positive integral powers
                            for (int e_index = 1; e_index < (int)lights[curr_light].pf; e_index++)
                                dpsl_exp *= dpsl;

                            // now dpsl_exp holds (dpsl)^pf power which is of course (s . l)^pf

                            i = 128 * dp * dpsl_exp / (atten);

                            r_sum2 += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                            g_sum2 += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                            b_sum2 += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
                        } // end if

                    } // end if

                    //Write_Error("\nexiting spotlight2, sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,   r_sum1, g_sum1, b_sum1,  r_sum2, g_sum2, b_sum2);

                } // end if spot light

            } // end for light

            // make sure colors aren't out of range
            if (r_sum0 > 255)
                r_sum0 = 255;
            if (g_sum0 > 255)
                g_sum0 = 255;
            if (b_sum0 > 255)
                b_sum0 = 255;

            if (r_sum1 > 255)
                r_sum1 = 255;
            if (g_sum1 > 255)
                g_sum1 = 255;
            if (b_sum1 > 255)
                b_sum1 = 255;

            if (r_sum2 > 255)
                r_sum2 = 255;
            if (g_sum2 > 255)
                g_sum2 = 255;
            if (b_sum2 > 255)
                b_sum2 = 255;

            //Write_Error("\nwriting color for poly %d", poly);

            //Write_Error("\n******** final sums rgb0[%d, %d, %d], rgb1[%d,%d,%d], rgb2[%d,%d,%d]", r_sum0, g_sum0, b_sum0,  r_sum1, g_sum1, b_sum1, r_sum2, g_sum2, b_sum2);

            // write the colors
            curr_poly->lit_color[0] = RGB16Bit565(r_sum0, g_sum0, b_sum0);
            curr_poly->lit_color[1] = RGB16Bit565(r_sum1, g_sum1, b_sum1);
            curr_poly->lit_color[2] = RGB16Bit565(r_sum2, g_sum2, b_sum2);

        }    // end if
        else // assume POLY4DV2_ATTR_SHADE_MODE_CONSTANT
        {
            // emmisive shading only, do nothing
            // ...
            curr_poly->lit_color[0] = curr_poly->color;

            //Write_Error("\nentering constant shader, and exiting...");

        } // end if

    } // end for poly

    // return success
    return (1);

} // end Light_RENDERLIST4DV2_World16

void World_To_Camera_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                    CAM4DV1_PTR cam)
{
    // NOTE: this is a matrix based function
    // this function transforms each polygon in the global render list
    // to camera coordinates based on the sent camera transform matrix
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list
    // the conversion of an object into polygons probably would have
    // happened after object culling, local transforms, local to world
    // and backface culling, so the minimum number of polygons from
    // each object are in the list, note that the function assumes
    // that at LEAST the local to world transform has been called
    // and the polygon data is in the transformed list tvlist of
    // the POLYF4DV1 object

    // transform each polygon in the render list into camera coordinates
    // assumes the render list has already been transformed to world
    // coordinates and the result is in tvlist[] of each polygon object

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // transform the vertex by the mcam matrix within the camera
            // it better be valid!
            POINT4D presult; // hold result of each transformation

            // transform point
            Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex].v, &cam->mcam, &presult);

            // store result back
            VECTOR4D_COPY(&curr_poly->tvlist[vertex].v, &presult);
        } // end for vertex

    } // end for poly

} // end World_To_Camera_RENDERLIST4DV2

void Sort_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list, int sort_method = SORT_POLYLIST_AVGZ)
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
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV2_PTR), Compare_AvgZ_POLYF4DV2);
    }
    break;

    case SORT_POLYLIST_NEARZ: // - sorts on closest z vertex of each poly
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV2_PTR), Compare_NearZ_POLYF4DV2);
    }
    break;

    case SORT_POLYLIST_FARZ: //  - sorts on farthest z vertex of each poly
    {
        qsort((void *)rend_list->poly_ptrs, rend_list->num_polys, sizeof(POLYF4DV2_PTR), Compare_FarZ_POLYF4DV2);
    }
    break;

    default:
        break;
    } // end switch

} // end Sort_RENDERLIST4DV2

void Camera_To_Perspective_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                          CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function transforms each polygon in the global render list
    // into perspective coordinates, based on the
    // sent camera object,
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list

    // transform each polygon in the render list into camera coordinates
    // assumes the render list has already been transformed to world
    // coordinates and the result is in tvlist[] of each polygon object

    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            float z = curr_poly->tvlist[vertex].z;

            // transform the vertex by the view parameters in the camera
            curr_poly->tvlist[vertex].x = cam->view_dist * curr_poly->tvlist[vertex].x / z;
            curr_poly->tvlist[vertex].y = cam->view_dist * curr_poly->tvlist[vertex].y * cam->aspect_ratio / z;
            // z = z, so no change

            // not that we are NOT dividing by the homogenous w coordinate since
            // we are not using a matrix operation for this version of the function

        } // end for vertex

    } // end for poly

} // end Camera_To_Perspective_RENDERLIST4DV2

void Perspective_To_Screen_RENDERLIST4DV2(RENDERLIST4DV2_PTR rend_list,
                                          CAM4DV1_PTR cam)
{
    // NOTE: this is not a matrix based function
    // this function transforms the perspective coordinates of the render
    // list into screen coordinates, based on the sent viewport in the camera
    // assuming that the viewplane coordinates were normalized
    // you would use this function instead of the object based function
    // if you decided earlier in the pipeline to turn each object into
    // a list of polygons and then add them to the global render list
    // you would only call this function if you previously performed
    // a normalized perspective transform

    // transform each polygon in the render list from perspective to screen
    // coordinates assumes the render list has already been transformed
    // to normalized perspective coordinates and the result is in tvlist[]
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        // acquire current polygon
        POLYF4DV2_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV2_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV2_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV2_STATE_BACKFACE))
            continue; // move onto next poly

        float alpha = (0.5 * cam->viewport_width - 0.5);
        float beta = (0.5 * cam->viewport_height - 0.5);

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // the vertex is in perspective normalized coords from -1 to 1
            // on each axis, simple scale them and invert y axis and project
            // to screen

            // transform the vertex by the view parameters in the camera
            curr_poly->tvlist[vertex].x = alpha + alpha * curr_poly->tvlist[vertex].x;
            curr_poly->tvlist[vertex].y = beta - beta * curr_poly->tvlist[vertex].y;
        } // end for vertex

    } // end for poly

} // end Perspective_To_Screen_RENDERLIST4DV2
int Compare_AvgZ_POLYF4DV2(const void *arg1, const void *arg2)
{
    // this function comapares the average z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV2_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV2_PTR *)(arg1));
    poly_2 = *((POLYF4DV2_PTR *)(arg2));

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

} // end Compare_AvgZ_POLYF4DV2

////////////////////////////////////////////////////////////////////////////////

int Compare_NearZ_POLYF4DV2(const void *arg1, const void *arg2)
{
    // this function comapares the closest z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV2_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV2_PTR *)(arg1));
    poly_2 = *((POLYF4DV2_PTR *)(arg2));

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

} // end Compare_NearZ_POLYF4DV2

////////////////////////////////////////////////////////////////////////////////

int Compare_FarZ_POLYF4DV2(const void *arg1, const void *arg2)
{
    // this function comapares the farthest z's of two polygons and is used by the
    // depth sort surface ordering algorithm

    float z1, z2;

    POLYF4DV2_PTR poly_1, poly_2;

    // dereference the poly pointers
    poly_1 = *((POLYF4DV2_PTR *)(arg1));
    poly_2 = *((POLYF4DV2_PTR *)(arg2));

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

} // end Compare_FarZ_POLYF4DV2

void DrawPhongTriangle(device_t *device, POLYF4DV2_PTR face) 
{
    int v0 = 0,
        v1 = 1,
        v2 = 2,
        temp = 0,
        tri_type = TRI_TYPE_NONE,
        irestart = INTERP_LHS;

    int dx, dy, dyl, dyr, // general deltas
        u, v, w,
        du, dv, dw,
        xi, yi,           // the current interpolated x,y
        ui, vi, wi,       // the current interpolated u,v
        index_x, index_y, // looping vars
        x, y,             // hold general x,y
        xstart,
        xend,
        ystart,
        yrestart,
        yend,
        xl,
        dxdyl,
        xr,
        dxdyr,
        dudyl,
        ul,
        dvdyl,
        vl,
        dwdyl,
        wl,
        dudyr,
        ur,
        dvdyr,
        vr,
        dwdyr,
        wr;

    int x0, y0, tu0, tv0, tw0, // cached vertices
        x1, y1, tu1, tv1, tw1,
        x2, y2, tu2, tv2, tw2;

    int r_base0, g_base0, b_base0,
        r_base1, g_base1, b_base1,
        r_base2, g_base2, b_base2;

    //USHORT *screen_ptr = NULL,
    //       *screen_line = NULL,
    //       *textmap = NULL,
    //       *dest_buffer = (USHORT *)_dest_buffer;

#ifdef DEBUG_ON
    // track rendering stats
    debug_polys_rendered_per_frame++;
#endif

    // adjust memory pitch to words, divide by 2

    // first trivial clipping rejection tests
    if (((face->tvlist[0].y < min_clip_y) &&
         (face->tvlist[1].y < min_clip_y) &&
         (face->tvlist[2].y < min_clip_y)) ||

        ((face->tvlist[0].y > max_clip_y) &&
         (face->tvlist[1].y > max_clip_y) &&
         (face->tvlist[2].y > max_clip_y)) ||

        ((face->tvlist[0].x < min_clip_x) &&
         (face->tvlist[1].x < min_clip_x) &&
         (face->tvlist[2].x < min_clip_x)) ||

        ((face->tvlist[0].x > max_clip_x) &&
         (face->tvlist[1].x > max_clip_x) &&
         (face->tvlist[2].x > max_clip_x)))
        return;

    // degenerate triangle
    if (((face->tvlist[0].x == face->tvlist[1].x) && (face->tvlist[1].x == face->tvlist[2].x)) ||
        ((face->tvlist[0].y == face->tvlist[1].y) && (face->tvlist[1].y == face->tvlist[2].y)))
        return;

    // v根据y从小到大排序， y0 < y1 < y2
    if (face->tvlist[v1].y < face->tvlist[v0].y)
    {
        SWAP(v0, v1, temp);
    }

    if (face->tvlist[v2].y < face->tvlist[v0].y)
    {
        SWAP(v0, v2, temp);
    }

    if (face->tvlist[v2].y < face->tvlist[v1].y)
    {
        SWAP(v1, v2, temp);
    }

    // now test for trivial flat sided cases
    if (face->tvlist[v0].y == face->tvlist[v1].y)
    {
        // set triangle type
        tri_type = TRI_TYPE_FLAT_TOP;

        // sort vertices left to right
        if (face->tvlist[v1].x < face->tvlist[v0].x)
        {
            SWAP(v0, v1, temp);
        }

    } // end if
    else
        // now test for trivial flat sided cases
        if (face->tvlist[v1].y == face->tvlist[v2].y)
    {
        // set triangle type
        tri_type = TRI_TYPE_FLAT_BOTTOM;

        // sort vertices left to right
        if (face->tvlist[v2].x < face->tvlist[v1].x)
        {
            SWAP(v1, v2, temp);
        }

    } // end if
    else
    {
        // must be a general triangle
        tri_type = TRI_TYPE_GENERAL;

    } // end else

    // assume 5.6.5 format -- sorry!
    // we can't afford a function call in the inner loops, so we must write
    // two hard coded versions, if we want support for both 5.6.5, and 5.5.5
    _RGB565FROM16BIT(face->lit_color[v0], &r_base0, &g_base0, &b_base0);
    _RGB565FROM16BIT(face->lit_color[v1], &r_base1, &g_base1, &b_base1);
    _RGB565FROM16BIT(face->lit_color[v2], &r_base2, &g_base2, &b_base2);

    // scale to 8 bit
    r_base0 <<= 3;
    g_base0 <<= 2;
    b_base0 <<= 3;

    // scale to 8 bit
    r_base1 <<= 3;
    g_base1 <<= 2;
    b_base1 <<= 3;

    // scale to 8 bit
    r_base2 <<= 3;
    g_base2 <<= 2;
    b_base2 <<= 3;

    // extract vertices for processing, now that we have order
    x0 = (int)(face->tvlist[v0].x + 0.5);
    y0 = (int)(face->tvlist[v0].y + 0.5);

    tu0 = r_base0;
    tv0 = g_base0;
    tw0 = b_base0;




    x1 = (int)(face->tvlist[v1].x + 0.5);
    y1 = (int)(face->tvlist[v1].y + 0.5);

    tu1 = r_base1;
    tv1 = g_base1;
    tw1 = b_base1;



    x2 = (int)(face->tvlist[v2].x + 0.5);
    y2 = (int)(face->tvlist[v2].y + 0.5);

    tu2 = r_base2;
    tv2 = g_base2;
    tw2 = b_base2;

    // set interpolation restart value
    yrestart = y1;

    // what kind of triangle
    if (tri_type & TRI_TYPE_FLAT_MASK)
    {

        if (tri_type == TRI_TYPE_FLAT_TOP)
        {
            // compute all deltas
            dy = (y2 - y0);

            dxdyl = ((x2 - x0) << FIXP16_SHIFT) / dy;
            dudyl = ((tu2 - tu0) << FIXP16_SHIFT) / dy;
            dvdyl = ((tv2 - tv0) << FIXP16_SHIFT) / dy;
            dwdyl = ((tw2 - tw0) << FIXP16_SHIFT) / dy;

            dxdyr = ((x2 - x1) << FIXP16_SHIFT) / dy;
            dudyr = ((tu2 - tu1) << FIXP16_SHIFT) / dy;
            dvdyr = ((tv2 - tv1) << FIXP16_SHIFT) / dy;
            dwdyr = ((tw2 - tw1) << FIXP16_SHIFT) / dy;

            // test for y clipping
            if (y0 < min_clip_y)
            {
                // compute overclip
                dy = (min_clip_y - y0);

                // computer new LHS starting values
                xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
                ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
                vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
                wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

                // compute new RHS starting values
                xr = dxdyr * dy + (x1 << FIXP16_SHIFT);
                ur = dudyr * dy + (tu1 << FIXP16_SHIFT);
                vr = dvdyr * dy + (tv1 << FIXP16_SHIFT);
                wr = dwdyr * dy + (tw1 << FIXP16_SHIFT);

                // compute new starting y
                ystart = min_clip_y;

            } // end if
            else
            {
                // no clipping

                // set starting values
                xl = (x0 << FIXP16_SHIFT);
                xr = (x1 << FIXP16_SHIFT);

                ul = (tu0 << FIXP16_SHIFT);
                vl = (tv0 << FIXP16_SHIFT);
                wl = (tw0 << FIXP16_SHIFT);

                ur = (tu1 << FIXP16_SHIFT);
                vr = (tv1 << FIXP16_SHIFT);
                wr = (tw1 << FIXP16_SHIFT);

                // set starting y
                ystart = y0;

            } // end else

        } // end if flat top
        else
        {
            // must be flat bottom

            // compute all deltas
            dy = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dy;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dy;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dy;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dy;

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dy;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dy;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dy;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dy;

            // test for y clipping
            if (y0 < min_clip_y)
            {
                // compute overclip
                dy = (min_clip_y - y0);

                // computer new LHS starting values
                xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
                ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
                vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
                wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

                // compute new RHS starting values
                xr = dxdyr * dy + (x0 << FIXP16_SHIFT);
                ur = dudyr * dy + (tu0 << FIXP16_SHIFT);
                vr = dvdyr * dy + (tv0 << FIXP16_SHIFT);
                wr = dwdyr * dy + (tw0 << FIXP16_SHIFT);

                // compute new starting y
                ystart = min_clip_y;

            } // end if
            else
            {
                // no clipping

                // set starting values
                xl = (x0 << FIXP16_SHIFT);
                xr = (x0 << FIXP16_SHIFT);

                ul = (tu0 << FIXP16_SHIFT);
                vl = (tv0 << FIXP16_SHIFT);
                wl = (tw0 << FIXP16_SHIFT);

                ur = (tu0 << FIXP16_SHIFT);
                vr = (tv0 << FIXP16_SHIFT);
                wr = (tw0 << FIXP16_SHIFT);

                // set starting y
                ystart = y0;

            } // end else

        } // end else flat bottom

        // test for bottom clip, always
        if ((yend = y2) > max_clip_y)
            yend = max_clip_y;

        // test for horizontal clipping
        if ((x0 < min_clip_x) || (x0 > max_clip_x) ||
            (x1 < min_clip_x) || (x1 > max_clip_x) ||
            (x2 < min_clip_x) || (x2 > max_clip_x))
        {
            // clip version

            std::cout<<"special clip"<<std::endl;
            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                ///////////////////////////////////////////////////////////////////////

                // test for x clipping, LHS
                if (xstart < min_clip_x)
                {
                    // compute x overlap
                    dx = min_clip_x - xstart;

                    // slide interpolants over
                    ui += dx * du;
                    vi += dx * dv;
                    wi += dx * dw;

                    // reset vars
                    xstart = min_clip_x;

                } // end if

                // test for x clipping RHS
                if (xend > max_clip_x)
                    xend = max_clip_x;

                ///////////////////////////////////////////////////////////////////////

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    // IUINT32 c = (ui<<16) | (vi<<8) | wi;
                    // IUINT32 c = (red << 16) | (green << 8) | blue;
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr

            } // end for y

        } // end if clip
        else
        {
            // non-clip version

            std::cout<<"special non-clip"<<std::endl;
            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi; 
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    // device_pixel(device, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

            } // end for y

        } // end if non-clipped

    } // end if
    else if (tri_type == TRI_TYPE_GENERAL)
    {

        // first test for bottom clip, always
        if ((yend = y2) > max_clip_y)
            yend = max_clip_y;

        // pre-test y clipping status
        if (y1 < min_clip_y)
        {
            // compute all deltas
            // LHS
            dyl = (y2 - y1);

            dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // compute overclip
            dyr = (min_clip_y - y0);
            dyl = (min_clip_y - y1);

            // computer new LHS starting values
            xl = dxdyl * dyl + (x1 << FIXP16_SHIFT);

            ul = dudyl * dyl + (tu1 << FIXP16_SHIFT);
            vl = dvdyl * dyl + (tv1 << FIXP16_SHIFT);
            wl = dwdyl * dyl + (tw1 << FIXP16_SHIFT);

            // compute new RHS starting values
            xr = dxdyr * dyr + (x0 << FIXP16_SHIFT);

            ur = dudyr * dyr + (tu0 << FIXP16_SHIFT);
            vr = dvdyr * dyr + (tv0 << FIXP16_SHIFT);
            wr = dwdyr * dyr + (tw0 << FIXP16_SHIFT);

            // compute new starting y
            ystart = min_clip_y;

            // test if we need swap to keep rendering left to right
            if (dxdyr > dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end if
        else if (y0 < min_clip_y)
        {
            // compute all deltas
            // LHS
            dyl = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // compute overclip
            dy = (min_clip_y - y0);

            // computer new LHS starting values
            xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
            ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
            vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
            wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

            // compute new RHS starting values
            xr = dxdyr * dy + (x0 << FIXP16_SHIFT);
            ur = dudyr * dy + (tu0 << FIXP16_SHIFT);
            vr = dvdyr * dy + (tv0 << FIXP16_SHIFT);
            wr = dwdyr * dy + (tw0 << FIXP16_SHIFT);

            // compute new starting y
            ystart = min_clip_y;

            // test if we need swap to keep rendering left to right
            if (dxdyr < dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end if
        else
        {
            // no initial y clipping

            // compute all deltas
            // LHS
            dyl = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // no clipping y

            // set starting values
            xl = (x0 << FIXP16_SHIFT);
            xr = (x0 << FIXP16_SHIFT);

            ul = (tu0 << FIXP16_SHIFT);
            vl = (tv0 << FIXP16_SHIFT);
            wl = (tw0 << FIXP16_SHIFT);

            ur = (tu0 << FIXP16_SHIFT);
            vr = (tv0 << FIXP16_SHIFT);
            wr = (tw0 << FIXP16_SHIFT);

            // set starting y
            ystart = y0;

            // test if we need swap to keep rendering left to right
            if (dxdyr < dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end else

        // test for horizontal clipping
        if ((x0 < min_clip_x) || (x0 > max_clip_x) ||
            (x1 < min_clip_x) || (x1 > max_clip_x) ||
            (x2 < min_clip_x) || (x2 > max_clip_x))
        {
            // clip version
            // x clipping

            std::cout<<"general clip"<<std::endl;
            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                ///////////////////////////////////////////////////////////////////////

                // test for x clipping, LHS
                if (xstart < min_clip_x)
                {
                    // compute x overlap
                    dx = min_clip_x - xstart;

                    // slide interpolants over
                    ui += dx * du;
                    vi += dx * dv;
                    wi += dx * dw;

                    // set x to left clip edge
                    xstart = min_clip_x;

                } // end if

                // test for x clipping RHS
                if (xend > max_clip_x)
                    xend = max_clip_x;

                ///////////////////////////////////////////////////////////////////////

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi; 
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    // device_pixel(device, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

                // test for yi hitting second region, if so change interpolant
                if (yi == yrestart)
                {
                    // test interpolation side change flag

                    if (irestart == INTERP_LHS)
                    {
                        // LHS
                        dyl = (y2 - y1);

                        dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
                        dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
                        dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
                        dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

                        // set starting values
                        xl = (x1 << FIXP16_SHIFT);
                        ul = (tu1 << FIXP16_SHIFT);
                        vl = (tv1 << FIXP16_SHIFT);
                        wl = (tw1 << FIXP16_SHIFT);

                        // interpolate down on LHS to even up
                        xl += dxdyl;
                        ul += dudyl;
                        vl += dvdyl;
                        wl += dwdyl;
                    } // end if
                    else
                    {
                        // RHS
                        dyr = (y1 - y2);

                        dxdyr = ((x1 - x2) << FIXP16_SHIFT) / dyr;
                        dudyr = ((tu1 - tu2) << FIXP16_SHIFT) / dyr;
                        dvdyr = ((tv1 - tv2) << FIXP16_SHIFT) / dyr;
                        dwdyr = ((tw1 - tw2) << FIXP16_SHIFT) / dyr;

                        // set starting values
                        xr = (x2 << FIXP16_SHIFT);
                        ur = (tu2 << FIXP16_SHIFT);
                        vr = (tv2 << FIXP16_SHIFT);
                        wr = (tw2 << FIXP16_SHIFT);

                        // interpolate down on RHS to even up
                        xr += dxdyr;
                        ur += dudyr;
                        vr += dvdyr;
                        wr += dwdyr;

                    } // end else

                } // end if

            } // end for y

        } // end if
        else
        {
            std::cout<<"general non-clip"<<std::endl;
            // no x clipping
            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi; 
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    // device_pixel(device, xi,  yi,  c);

                    //将颜色插值，改成发现和pos的插值
                    //得到normal和pos之后，算光照

                    VECTOR4D fragNormal;
                    VECTOR4D fragPos;


                    //得到poly原本的颜色
                    int r_base, g_base, b_base;
                    _RGB565FROM16BIT(face->color, &r_base, &g_base, &b_base);
                    r_base <<= 3;
                    g_base <<= 2;
                    b_base <<= 3;

                    int color;
                    CAM4DV1 cam;

                    //计算出该像素的最终光照颜色
                    ComputePhongShadingPixelColor(r_base, g_base, b_base, lights, &cam, &fragPos, &fragNormal, color);

                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);

                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

                // test for yi hitting second region, if so change interpolant
                // y到达了转折点，则更新斜率
                if (yi == yrestart)
                {
                    // test interpolation side change flag

                    if (irestart == INTERP_LHS)
                    {
                        // LHS
                        dyl = (y2 - y1);

						

                        dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
                        dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
                        dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
                        dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

                        // set starting values
                        xl = (x1 << FIXP16_SHIFT);
                        ul = (tu1 << FIXP16_SHIFT);
                        vl = (tv1 << FIXP16_SHIFT);
                        wl = (tw1 << FIXP16_SHIFT);

                        // interpolate down on LHS to even up
                        xl += dxdyl;
                        ul += dudyl;
                        vl += dvdyl;
                        wl += dwdyl;
                    } // end if
                    else
                    {
                        // RHS

                        dyr = (y1 - y2);

                        dxdyr = ((x1 - x2) << FIXP16_SHIFT) / dyr;
                        dudyr = ((tu1 - tu2) << FIXP16_SHIFT) / dyr;
                        dvdyr = ((tv1 - tv2) << FIXP16_SHIFT) / dyr;
                        dwdyr = ((tw1 - tw2) << FIXP16_SHIFT) / dyr;

                        // set starting values
                        xr = (x2 << FIXP16_SHIFT);
                        ur = (tu2 << FIXP16_SHIFT);
                        vr = (tv2 << FIXP16_SHIFT);
                        wr = (tw2 << FIXP16_SHIFT);

                        // interpolate down on RHS to even up
                        xr += dxdyr;
                        ur += dudyr;
                        vr += dvdyr;
                        wr += dwdyr;
                    } // end else

                } // end if

            } // end for y

        } // end else

    } // end if


}

void Draw_Gouraud_Triangle16(device_t *device, POLYF4DV2_PTR face) 
{
    // this function draws a gouraud shaded polygon, based on the affine texture mapper, instead
    // of interpolating the texture coordinates, we simply interpolate the (R,G,B) values across
    // the polygons, I simply needed at another interpolant, I have mapped u->red, v->green, w->blue

    int v0 = 0,
        v1 = 1,
        v2 = 2,
        temp = 0,
        tri_type = TRI_TYPE_NONE,
        irestart = INTERP_LHS;

    int dx, dy, dyl, dyr, // general deltas
        u, v, w,
        du, dv, dw,
        xi, yi,           // the current interpolated x,y
        ui, vi, wi,       // the current interpolated u,v
        index_x, index_y, // looping vars
        x, y,             // hold general x,y
        xstart,
        xend,
        ystart,
        yrestart,
        yend,
        xl,
        dxdyl,
        xr,
        dxdyr,
        dudyl,
        ul,
        dvdyl,
        vl,
        dwdyl,
        wl,
        dudyr,
        ur,
        dvdyr,
        vr,
        dwdyr,
        wr;

    int x0, y0, tu0, tv0, tw0, // cached vertices
        x1, y1, tu1, tv1, tw1,
        x2, y2, tu2, tv2, tw2;

    int r_base0, g_base0, b_base0,
        r_base1, g_base1, b_base1,
        r_base2, g_base2, b_base2;

    //USHORT *screen_ptr = NULL,
    //       *screen_line = NULL,
    //       *textmap = NULL,
    //       *dest_buffer = (USHORT *)_dest_buffer;

#ifdef DEBUG_ON
    // track rendering stats
    debug_polys_rendered_per_frame++;
#endif

    // adjust memory pitch to words, divide by 2

    // first trivial clipping rejection tests
    if (((face->tvlist[0].y < min_clip_y) &&
         (face->tvlist[1].y < min_clip_y) &&
         (face->tvlist[2].y < min_clip_y)) ||

        ((face->tvlist[0].y > max_clip_y) &&
         (face->tvlist[1].y > max_clip_y) &&
         (face->tvlist[2].y > max_clip_y)) ||

        ((face->tvlist[0].x < min_clip_x) &&
         (face->tvlist[1].x < min_clip_x) &&
         (face->tvlist[2].x < min_clip_x)) ||

        ((face->tvlist[0].x > max_clip_x) &&
         (face->tvlist[1].x > max_clip_x) &&
         (face->tvlist[2].x > max_clip_x)))
        return;

    // degenerate triangle
    if (((face->tvlist[0].x == face->tvlist[1].x) && (face->tvlist[1].x == face->tvlist[2].x)) ||
        ((face->tvlist[0].y == face->tvlist[1].y) && (face->tvlist[1].y == face->tvlist[2].y)))
        return;

    // sort vertices
    if (face->tvlist[v1].y < face->tvlist[v0].y)
    {
        SWAP(v0, v1, temp);
    }

    if (face->tvlist[v2].y < face->tvlist[v0].y)
    {
        SWAP(v0, v2, temp);
    }

    if (face->tvlist[v2].y < face->tvlist[v1].y)
    {
        SWAP(v1, v2, temp);
    }

    // now test for trivial flat sided cases
    if (face->tvlist[v0].y == face->tvlist[v1].y)
    {
        // set triangle type
        tri_type = TRI_TYPE_FLAT_TOP;

        // sort vertices left to right
        if (face->tvlist[v1].x < face->tvlist[v0].x)
        {
            SWAP(v0, v1, temp);
        }

    } // end if
    else
        // now test for trivial flat sided cases
        if (face->tvlist[v1].y == face->tvlist[v2].y)
    {
        // set triangle type
        tri_type = TRI_TYPE_FLAT_BOTTOM;

        // sort vertices left to right
        if (face->tvlist[v2].x < face->tvlist[v1].x)
        {
            SWAP(v1, v2, temp);
        }

    } // end if
    else
    {
        // must be a general triangle
        tri_type = TRI_TYPE_GENERAL;

    } // end else

    // assume 5.6.5 format -- sorry!
    // we can't afford a function call in the inner loops, so we must write
    // two hard coded versions, if we want support for both 5.6.5, and 5.5.5
    _RGB565FROM16BIT(face->lit_color[v0], &r_base0, &g_base0, &b_base0);
    _RGB565FROM16BIT(face->lit_color[v1], &r_base1, &g_base1, &b_base1);
    _RGB565FROM16BIT(face->lit_color[v2], &r_base2, &g_base2, &b_base2);

    // scale to 8 bit
    r_base0 <<= 3;
    g_base0 <<= 2;
    b_base0 <<= 3;

    // scale to 8 bit
    r_base1 <<= 3;
    g_base1 <<= 2;
    b_base1 <<= 3;

    // scale to 8 bit
    r_base2 <<= 3;
    g_base2 <<= 2;
    b_base2 <<= 3;

    // extract vertices for processing, now that we have order
    x0 = (int)(face->tvlist[v0].x + 0.5);
    y0 = (int)(face->tvlist[v0].y + 0.5);

    tu0 = r_base0;
    tv0 = g_base0;
    tw0 = b_base0;

    x1 = (int)(face->tvlist[v1].x + 0.5);
    y1 = (int)(face->tvlist[v1].y + 0.5);

    tu1 = r_base1;
    tv1 = g_base1;
    tw1 = b_base1;

    x2 = (int)(face->tvlist[v2].x + 0.5);
    y2 = (int)(face->tvlist[v2].y + 0.5);

    tu2 = r_base2;
    tv2 = g_base2;
    tw2 = b_base2;

    // set interpolation restart value
    yrestart = y1;

    // what kind of triangle
    if (tri_type & TRI_TYPE_FLAT_MASK)
    {

        if (tri_type == TRI_TYPE_FLAT_TOP)
        {
            // compute all deltas
            dy = (y2 - y0);

            dxdyl = ((x2 - x0) << FIXP16_SHIFT) / dy;
            dudyl = ((tu2 - tu0) << FIXP16_SHIFT) / dy;
            dvdyl = ((tv2 - tv0) << FIXP16_SHIFT) / dy;
            dwdyl = ((tw2 - tw0) << FIXP16_SHIFT) / dy;

            dxdyr = ((x2 - x1) << FIXP16_SHIFT) / dy;
            dudyr = ((tu2 - tu1) << FIXP16_SHIFT) / dy;
            dvdyr = ((tv2 - tv1) << FIXP16_SHIFT) / dy;
            dwdyr = ((tw2 - tw1) << FIXP16_SHIFT) / dy;

            // test for y clipping
            if (y0 < min_clip_y)
            {
                // compute overclip
                dy = (min_clip_y - y0);

                // computer new LHS starting values
                xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
                ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
                vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
                wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

                // compute new RHS starting values
                xr = dxdyr * dy + (x1 << FIXP16_SHIFT);
                ur = dudyr * dy + (tu1 << FIXP16_SHIFT);
                vr = dvdyr * dy + (tv1 << FIXP16_SHIFT);
                wr = dwdyr * dy + (tw1 << FIXP16_SHIFT);

                // compute new starting y
                ystart = min_clip_y;

            } // end if
            else
            {
                // no clipping

                // set starting values
                xl = (x0 << FIXP16_SHIFT);
                xr = (x1 << FIXP16_SHIFT);

                ul = (tu0 << FIXP16_SHIFT);
                vl = (tv0 << FIXP16_SHIFT);
                wl = (tw0 << FIXP16_SHIFT);

                ur = (tu1 << FIXP16_SHIFT);
                vr = (tv1 << FIXP16_SHIFT);
                wr = (tw1 << FIXP16_SHIFT);

                // set starting y
                ystart = y0;

            } // end else

        } // end if flat top
        else
        {
            // must be flat bottom

            // compute all deltas
            dy = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dy;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dy;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dy;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dy;

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dy;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dy;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dy;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dy;

            // test for y clipping
            if (y0 < min_clip_y)
            {
                // compute overclip
                dy = (min_clip_y - y0);

                // computer new LHS starting values
                xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
                ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
                vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
                wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

                // compute new RHS starting values
                xr = dxdyr * dy + (x0 << FIXP16_SHIFT);
                ur = dudyr * dy + (tu0 << FIXP16_SHIFT);
                vr = dvdyr * dy + (tv0 << FIXP16_SHIFT);
                wr = dwdyr * dy + (tw0 << FIXP16_SHIFT);

                // compute new starting y
                ystart = min_clip_y;

            } // end if
            else
            {
                // no clipping

                // set starting values
                xl = (x0 << FIXP16_SHIFT);
                xr = (x0 << FIXP16_SHIFT);

                ul = (tu0 << FIXP16_SHIFT);
                vl = (tv0 << FIXP16_SHIFT);
                wl = (tw0 << FIXP16_SHIFT);

                ur = (tu0 << FIXP16_SHIFT);
                vr = (tv0 << FIXP16_SHIFT);
                wr = (tw0 << FIXP16_SHIFT);

                // set starting y
                ystart = y0;

            } // end else

        } // end else flat bottom

        // test for bottom clip, always
        if ((yend = y2) > max_clip_y)
            yend = max_clip_y;

        // test for horizontal clipping
        if ((x0 < min_clip_x) || (x0 > max_clip_x) ||
            (x1 < min_clip_x) || (x1 > max_clip_x) ||
            (x2 < min_clip_x) || (x2 > max_clip_x))
        {
            // clip version

            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                ///////////////////////////////////////////////////////////////////////

                // test for x clipping, LHS
                if (xstart < min_clip_x)
                {
                    // compute x overlap
                    dx = min_clip_x - xstart;

                    // slide interpolants over
                    ui += dx * du;
                    vi += dx * dv;
                    wi += dx * dw;

                    // reset vars
                    xstart = min_clip_x;

                } // end if

                // test for x clipping RHS
                if (xend > max_clip_x)
                    xend = max_clip_x;

                ///////////////////////////////////////////////////////////////////////

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    // IUINT32 c = (ui<<16) | (vi<<8) | wi;
                    // IUINT32 c = (red << 16) | (green << 8) | blue;
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr

            } // end for y

        } // end if clip
        else
        {
            // non-clip version

            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi; 
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    // device_pixel(device, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

            } // end for y

        } // end if non-clipped

    } // end if
    else if (tri_type == TRI_TYPE_GENERAL)
    {

        // first test for bottom clip, always
        if ((yend = y2) > max_clip_y)
            yend = max_clip_y;

        // pre-test y clipping status
        if (y1 < min_clip_y)
        {
            // compute all deltas
            // LHS
            dyl = (y2 - y1);

            dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // compute overclip
            dyr = (min_clip_y - y0);
            dyl = (min_clip_y - y1);

            // computer new LHS starting values
            xl = dxdyl * dyl + (x1 << FIXP16_SHIFT);

            ul = dudyl * dyl + (tu1 << FIXP16_SHIFT);
            vl = dvdyl * dyl + (tv1 << FIXP16_SHIFT);
            wl = dwdyl * dyl + (tw1 << FIXP16_SHIFT);

            // compute new RHS starting values
            xr = dxdyr * dyr + (x0 << FIXP16_SHIFT);

            ur = dudyr * dyr + (tu0 << FIXP16_SHIFT);
            vr = dvdyr * dyr + (tv0 << FIXP16_SHIFT);
            wr = dwdyr * dyr + (tw0 << FIXP16_SHIFT);

            // compute new starting y
            ystart = min_clip_y;

            // test if we need swap to keep rendering left to right
            if (dxdyr > dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end if
        else if (y0 < min_clip_y)
        {
            // compute all deltas
            // LHS
            dyl = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // compute overclip
            dy = (min_clip_y - y0);

            // computer new LHS starting values
            xl = dxdyl * dy + (x0 << FIXP16_SHIFT);
            ul = dudyl * dy + (tu0 << FIXP16_SHIFT);
            vl = dvdyl * dy + (tv0 << FIXP16_SHIFT);
            wl = dwdyl * dy + (tw0 << FIXP16_SHIFT);

            // compute new RHS starting values
            xr = dxdyr * dy + (x0 << FIXP16_SHIFT);
            ur = dudyr * dy + (tu0 << FIXP16_SHIFT);
            vr = dvdyr * dy + (tv0 << FIXP16_SHIFT);
            wr = dwdyr * dy + (tw0 << FIXP16_SHIFT);

            // compute new starting y
            ystart = min_clip_y;

            // test if we need swap to keep rendering left to right
            if (dxdyr < dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end if
        else
        {
            // no initial y clipping

            // compute all deltas
            // LHS
            dyl = (y1 - y0);

            dxdyl = ((x1 - x0) << FIXP16_SHIFT) / dyl;
            dudyl = ((tu1 - tu0) << FIXP16_SHIFT) / dyl;
            dvdyl = ((tv1 - tv0) << FIXP16_SHIFT) / dyl;
            dwdyl = ((tw1 - tw0) << FIXP16_SHIFT) / dyl;

            // RHS
            dyr = (y2 - y0);

            dxdyr = ((x2 - x0) << FIXP16_SHIFT) / dyr;
            dudyr = ((tu2 - tu0) << FIXP16_SHIFT) / dyr;
            dvdyr = ((tv2 - tv0) << FIXP16_SHIFT) / dyr;
            dwdyr = ((tw2 - tw0) << FIXP16_SHIFT) / dyr;

            // no clipping y

            // set starting values
            xl = (x0 << FIXP16_SHIFT);
            xr = (x0 << FIXP16_SHIFT);

            ul = (tu0 << FIXP16_SHIFT);
            vl = (tv0 << FIXP16_SHIFT);
            wl = (tw0 << FIXP16_SHIFT);

            ur = (tu0 << FIXP16_SHIFT);
            vr = (tv0 << FIXP16_SHIFT);
            wr = (tw0 << FIXP16_SHIFT);

            // set starting y
            ystart = y0;

            // test if we need swap to keep rendering left to right
            if (dxdyr < dxdyl)
            {
                SWAP(dxdyl, dxdyr, temp);
                SWAP(dudyl, dudyr, temp);
                SWAP(dvdyl, dvdyr, temp);
                SWAP(dwdyl, dwdyr, temp);
                SWAP(xl, xr, temp);
                SWAP(ul, ur, temp);
                SWAP(vl, vr, temp);
                SWAP(wl, wr, temp);
                SWAP(x1, x2, temp);
                SWAP(y1, y2, temp);
                SWAP(tu1, tu2, temp);
                SWAP(tv1, tv2, temp);
                SWAP(tw1, tw2, temp);

                // set interpolation restart
                irestart = INTERP_RHS;

            } // end if

        } // end else

        // test for horizontal clipping
        if ((x0 < min_clip_x) || (x0 > max_clip_x) ||
            (x1 < min_clip_x) || (x1 > max_clip_x) ||
            (x2 < min_clip_x) || (x2 > max_clip_x))
        {
            // clip version
            // x clipping

            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                ///////////////////////////////////////////////////////////////////////

                // test for x clipping, LHS
                if (xstart < min_clip_x)
                {
                    // compute x overlap
                    dx = min_clip_x - xstart;

                    // slide interpolants over
                    ui += dx * du;
                    vi += dx * dv;
                    wi += dx * dw;

                    // set x to left clip edge
                    xstart = min_clip_x;

                } // end if

                // test for x clipping RHS
                if (xend > max_clip_x)
                    xend = max_clip_x;

                ///////////////////////////////////////////////////////////////////////

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi; 
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    // device_pixel(device, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

                // test for yi hitting second region, if so change interpolant
                if (yi == yrestart)
                {
                    // test interpolation side change flag

                    if (irestart == INTERP_LHS)
                    {
                        // LHS
                        dyl = (y2 - y1);

                        dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
                        dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
                        dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
                        dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

                        // set starting values
                        xl = (x1 << FIXP16_SHIFT);
                        ul = (tu1 << FIXP16_SHIFT);
                        vl = (tv1 << FIXP16_SHIFT);
                        wl = (tw1 << FIXP16_SHIFT);

                        // interpolate down on LHS to even up
                        xl += dxdyl;
                        ul += dudyl;
                        vl += dvdyl;
                        wl += dwdyl;
                    } // end if
                    else
                    {
                        // RHS
                        dyr = (y1 - y2);

                        dxdyr = ((x1 - x2) << FIXP16_SHIFT) / dyr;
                        dudyr = ((tu1 - tu2) << FIXP16_SHIFT) / dyr;
                        dvdyr = ((tv1 - tv2) << FIXP16_SHIFT) / dyr;
                        dwdyr = ((tw1 - tw2) << FIXP16_SHIFT) / dyr;

                        // set starting values
                        xr = (x2 << FIXP16_SHIFT);
                        ur = (tu2 << FIXP16_SHIFT);
                        vr = (tv2 << FIXP16_SHIFT);
                        wr = (tw2 << FIXP16_SHIFT);

                        // interpolate down on RHS to even up
                        xr += dxdyr;
                        ur += dudyr;
                        vr += dvdyr;
                        wr += dwdyr;

                    } // end else

                } // end if

            } // end for y

        } // end if
        else
        {
            // no x clipping
            // point screen ptr to starting line
            // screen_ptr = dest_buffer + (ystart * mem_pitch);

            for (yi = ystart; yi <= yend; yi++)
            {
                // compute span endpoints
                xstart = ((xl + FIXP16_ROUND_UP) >> FIXP16_SHIFT);
                xend = ((xr + FIXP16_ROUND_UP) >> FIXP16_SHIFT);

                // compute starting points for u,v,w interpolants
                ui = ul + FIXP16_ROUND_UP;
                vi = vl + FIXP16_ROUND_UP;
                wi = wl + FIXP16_ROUND_UP;

                // compute u,v interpolants
                if ((dx = (xend - xstart)) > 0)
                {
                    du = (ur - ul) / dx;
                    dv = (vr - vl) / dx;
                    dw = (wr - wl) / dx;
                } // end if
                else
                {
                    du = (ur - ul);
                    dv = (vr - vl);
                    dw = (wr - wl);
                } // end else

                // draw span
                for (xi = xstart; xi <= xend; xi++)
                {
                    // write textel assume 5.6.5
                    // screen_ptr[xi] = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));

                    // IUINT32 c = (ui<<16) | (vi<<8) | wi; 
                    // IUINT32 c = ((ui >> (FIXP16_SHIFT + 3))<<16) | ((vi >> (FIXP16_SHIFT + 2))<<8) | (wi >> (FIXP16_SHIFT + 3)); 
                    // device_pixel(device, xi,  yi,  c);

                    int color = ((ui >> (FIXP16_SHIFT + 3)) << 11) + ((vi >> (FIXP16_SHIFT + 2)) << 5) + (wi >> (FIXP16_SHIFT + 3));
                    IUINT32 c;
                    RGBFrom565(color, c);
                    device_pixel(device, xi,  yi,  c);
                    // interpolate u,v
                    ui += du;
                    vi += dv;
                    wi += dw;
                } // end for xi

                // interpolate u,v,w,x along right and left edge
                xl += dxdyl;
                ul += dudyl;
                vl += dvdyl;
                wl += dwdyl;

                xr += dxdyr;
                ur += dudyr;
                vr += dvdyr;
                wr += dwdyr;

                // advance screen ptr
                // screen_ptr += mem_pitch;

                // test for yi hitting second region, if so change interpolant
                if (yi == yrestart)
                {
                    // test interpolation side change flag

                    if (irestart == INTERP_LHS)
                    {
                        // LHS
                        dyl = (y2 - y1);

                        dxdyl = ((x2 - x1) << FIXP16_SHIFT) / dyl;
                        dudyl = ((tu2 - tu1) << FIXP16_SHIFT) / dyl;
                        dvdyl = ((tv2 - tv1) << FIXP16_SHIFT) / dyl;
                        dwdyl = ((tw2 - tw1) << FIXP16_SHIFT) / dyl;

                        // set starting values
                        xl = (x1 << FIXP16_SHIFT);
                        ul = (tu1 << FIXP16_SHIFT);
                        vl = (tv1 << FIXP16_SHIFT);
                        wl = (tw1 << FIXP16_SHIFT);

                        // interpolate down on LHS to even up
                        xl += dxdyl;
                        ul += dudyl;
                        vl += dvdyl;
                        wl += dwdyl;
                    } // end if
                    else
                    {
                        // RHS
                        dyr = (y1 - y2);

                        dxdyr = ((x1 - x2) << FIXP16_SHIFT) / dyr;
                        dudyr = ((tu1 - tu2) << FIXP16_SHIFT) / dyr;
                        dvdyr = ((tv1 - tv2) << FIXP16_SHIFT) / dyr;
                        dwdyr = ((tw1 - tw2) << FIXP16_SHIFT) / dyr;

                        // set starting values
                        xr = (x2 << FIXP16_SHIFT);
                        ur = (tu2 << FIXP16_SHIFT);
                        vr = (tv2 << FIXP16_SHIFT);
                        wr = (tw2 << FIXP16_SHIFT);

                        // interpolate down on RHS to even up
                        xr += dxdyr;
                        ur += dudyr;
                        vr += dvdyr;
                        wr += dwdyr;
                    } // end else

                } // end if

            } // end for y

        } // end else


    } // end if

} // end Draw_Gouraud_Triangle16

void ComputePhongShadingPixelColor(int r_base, int g_base, int b_base, LIGHTV1_PTR prtlights, CAM4DV1_PTR ptrCam, VECTOR4D_PTR prtFragPos, VECTOR4D_PTR ptrFragNormal, int &color)
{
    VECTOR4D lightPos;
    VECTOR4D lightDir; //光线方向
    LIGHTV1 ptrCurLight;
    int sumR, sumG, sumB;
    int max_lights = 4;

    for (int curr_light = 0; curr_light < max_lights; curr_light++)
    {
        if (lights[curr_light].state == LIGHTV1_STATE_OFF)
            continue;

        if (lights[curr_light].attr & LIGHTV1_ATTR_AMBIENT) //环境光
        {
            int ri,gi,bi;
            ri = ((lights[curr_light].c_ambient.r * r_base) / 256);
            gi = ((lights[curr_light].c_ambient.g * g_base) / 256);
            bi = ((lights[curr_light].c_ambient.b * b_base) / 256);

            sumR += ri;
            sumG += gi;
            sumB += bi;
        } // end if
        else if (lights[curr_light].attr & LIGHTV1_ATTR_INFINITE) //平行光
        {
            float dp = VECTOR4D_Dot(ptrFragNormal, &lights[curr_light].dir);
            if (dp > 0)
            {
                float i = 128 * dp;
                sumR += (lights[curr_light].c_diffuse.r * r_base * i) / (256 * 128);
                sumG += (lights[curr_light].c_diffuse.g * g_base * i) / (256 * 128);
                sumB += (lights[curr_light].c_diffuse.b * b_base * i) / (256 * 128);
            } // end if
        }
        else if (lights[curr_light].attr & LIGHTV1_ATTR_POINT) //点光源
        {
            VECTOR4D l; //光线方向
            VECTOR4D_Build(prtFragPos, &lights[curr_light].pos, &l);
            float dp = VECTOR4D_Dot(ptrFragNormal, &l);
            if (dp > 0)
            {
                sumR += (lights[curr_light].c_diffuse.r * dp * r_base)/256;
                sumG += (lights[curr_light].c_diffuse.g * dp * g_base)/256;
                sumB += (lights[curr_light].c_diffuse.b * dp * b_base)/256;
            } // end if
        }
    }

    sumR = max(255, sumR);
    sumG = max(255, sumG);
    sumB = max(255, sumB);

    color = ((sumR >> (FIXP16_SHIFT + 3)) << 11) + ((sumG >> (FIXP16_SHIFT + 2)) << 5) + (sumB >> (FIXP16_SHIFT + 3));
}