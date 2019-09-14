#include "frame.h"

// storage for our lookup tables
float cos_look[361]; // 1 extra element so we can store 0-360 inclusive
float sin_look[361];

USHORT RGB16Bit565(int r, int g, int b)
{
    // this function simply builds a 5.6.5 format 16 bit pixel
    // assumes input is RGB 0-255 each channel
    r >>= 3;
    g >>= 2;
    b >>= 3;
    return (_RGB16BIT565((r), (g), (b)));

} // end RGB16Bit565

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList)
{
    renderList->num_polys = 0;
}

void Init_CAM4DV1(CAM4DV1_PTR cam,        // the camera object
                  int cam_attr,           // attributes
                  POINT4D_PTR cam_pos,    // initial camera position
                  VECTOR4D_PTR cam_dir,   // initial camera angles
                  POINT4D_PTR cam_target, // UVN target
                  float near_clip_z,      // near and far clipping planes
                  float far_clip_z,
                  float fov,            // field of view in degrees
                  float viewport_width, // size of final screen viewport
                  float viewport_height)
{
    // this function initializes the camera object cam, the function
    // doesn't do a lot of error checking or sanity checking since
    // I want to allow you to create projections as you wish, also
    // I tried to minimize the number of parameters the functions needs

    // first set up parms that are no brainers
    cam->attr = cam_attr; // camera attributes

    VECTOR4D_COPY(&cam->pos, cam_pos); // positions
    VECTOR4D_COPY(&cam->dir, cam_dir); // direction vector or angles for
                                       // euler camera
    // for UVN camera
    VECTOR4D_INITXYZ(&cam->u, 1, 0, 0); // set to +x
    VECTOR4D_INITXYZ(&cam->v, 0, 1, 0); // set to +y
    VECTOR4D_INITXYZ(&cam->n, 0, 0, 1); // set to +z

    if (cam_target != nullptr)
        VECTOR4D_COPY(&cam->target, cam_target); // UVN target
    else
        VECTOR4D_ZERO(&cam->target);

    cam->near_clip_z = near_clip_z; // near z=constant clipping plane
    cam->far_clip_z = far_clip_z;   // far z=constant clipping plane

    cam->viewport_width = viewport_width; // dimensions of viewport
    cam->viewport_height = viewport_height;

    cam->viewport_center_x = (viewport_width - 1) / 2; // center of viewport
    cam->viewport_center_y = (viewport_height - 1) / 2;

    cam->aspect_ratio = (float)viewport_width / (float)viewport_height;

    // set all camera matrices to identity matrix
    MAT_IDENTITY_4X4(&cam->mcam);
    MAT_IDENTITY_4X4(&cam->mper);
    MAT_IDENTITY_4X4(&cam->mscr);

    // set independent vars
    cam->fov = fov;

    // set the viewplane dimensions up, they will be 2 x (2/ar)
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

void MAT_IDENTITY_4X4(MATRIX4X4_PTR m)
{
    memcpy((void *)(m), (void *)&IMAT_4X4, sizeof(MATRIX4X4));
}

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

float VECTOR3D_Length(VECTOR3D_PTR va)
{
    // computes the magnitude of a vector, slow

    return ((float)sqrtf(va->x * va->x + va->y * va->y + va->z * va->z));

} // end VECTOR3D_Length

int Insert_POLYF4DV1_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                    POLYF4DV1_PTR poly)
{
    // inserts the sent polyface POLYF4DV1 into the render list

    // step 0: are we full?
    if (rend_list->num_polys >= RENDERLIST4DV1_MAX_POLYS)
        return (0);

    // step 1: copy polygon into next opening in polygon render list

    // point pointer to polygon structure
    rend_list->poly_ptrs[rend_list->num_polys] = &rend_list->poly_data[rend_list->num_polys];

    // copy face right into array, thats it
    memcpy((void *)&rend_list->poly_data[rend_list->num_polys], (void *)poly, sizeof(POLYF4DV1));

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

} // end Insert_POLYF4DV1_RENDERLIST4DV1

void Build_XYZ_Rotation_MATRIX4X4(float theta_x, // euler angles
                                  float theta_y,
                                  float theta_z,
                                  MATRIX4X4_PTR mrot) // output
{
    // this helper function takes a set if euler angles and computes
    // a rotation matrix from them, usefull for object and camera
    // work, also  we will do a little testing in the function to determine
    // the rotations that need to be performed, since there's no
    // reason to perform extra matrix multiplies if the angles are
    // zero!

    MATRIX4X4 mx, my, mz, mtmp;         // working matrices
    float sin_theta = 0, cos_theta = 0; // used to initialize matrices
    int rot_seq = 0;                    // 1 for x, 2 for y, 4 for z

    // step 0: fill in with identity matrix
    MAT_IDENTITY_4X4(mrot);

    // step 1: based on zero and non-zero rotation angles, determine
    // rotation sequence
    if (fabs(theta_x) > EPSILON_E5) // x
        rot_seq = rot_seq | 1;

    if (fabs(theta_y) > EPSILON_E5) // y
        rot_seq = rot_seq | 2;

    if (fabs(theta_z) > EPSILON_E5) // z
        rot_seq = rot_seq | 4;

    // now case on sequence
    switch (rot_seq)
    {
    case 0: // no rotation
    {
        // what a waste!
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

        // that's it, copy to output matrix
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
        Mat_Init_4X4(&my, cos_theta, 0, -sin_theta, 0,
                     0, 1, 0, 0,
                     sin_theta, 0, cos_theta, 0,
                     0, 0, 0, 1);

        // that's it, copy to output matrix
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
    return (sin_look[theta_int] +
            theta_frac * (sin_look[theta_int + 1] - sin_look[theta_int]));

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
    return (cos_look[theta_int] +
            theta_frac * (cos_look[theta_int + 1] - cos_look[theta_int]));

} // end Fast_Cos

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

void Mat_Mul_4X4(MATRIX4X4_PTR ma,
                 MATRIX4X4_PTR mb,
                 MATRIX4X4_PTR mprod)
{
    // this function multiplies two 4x4 matrices together and
    // and stores the result in mprod
    // note later we will take advantage of the fact that we know
    // that w=1 always, and that the last column of a 4x4 is
    // always 0

    for (int row = 0; row < 4; row++)
    {
        for (int col = 0; col < 4; col++)
        {
            // compute dot product from row of ma
            // and column of mb

            float sum = 0; // used to hold result

            for (int index = 0; index < 4; index++)
            {
                // add in next product pair
                sum += (ma->M[row][index] * mb->M[index][col]);
            } // end for index

            // insert resulting row,col element
            mprod->M[row][col] = sum;

        } // end for col

    } // end for row

} // end Mat_Mul_4X4

void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, // render list to transform
                              MATRIX4X4_PTR mt,             // transformation matrix
                              int coord_select)             // selects coords to transform
{
    // this function simply transforms all of the polygons vertices in the local or trans
    // array of the render list by the sent matrix

    // what coordinates should be transformed?
    switch (coord_select)
    {
    case TRANSFORM_LOCAL_ONLY:
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) || (curr_poly->state & POLY4DV1_STATE_CLIPPED) || (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue;

            for (int vertex = 0; vertex < 3; vertex++)
            {
                POINT4D presult;

                // transform point
                Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &presult);

                VECTOR4D_COPY(&curr_poly->vlist[vertex], &presult);
            } // end for vertex

        } // end for poly
    }
    break;

    case TRANSFORM_TRANS_ONLY:
    {
        // transform each "transformed" vertex of the render list
        // remember, the idea of the tvlist[] array is to accumulate
        // transformations
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // is this polygon valid?
            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // transform the vertex by mt
                POINT4D presult; // hold result of each transformation

                // transform point
                Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], mt, &presult);

                // store result back
                VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
            } // end for vertex

        } // end for poly
    }
    break;

    case TRANSFORM_LOCAL_TO_TRANS:
    {
        // transform each local/model vertex of the render list and store result
        // in "transformed" vertex list
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // is this polygon valid?
            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // transform the vertex by mt
                Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &curr_poly->tvlist[vertex]);
            } // end for vertex

        } // end for poly
    }
    break;

    default:
        break;

    } // end switch

} // end Transform_RENDERLIST4DV1

void Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                   POINT4D_PTR world_pos,
                                   int coord_select)
{
    // NOTE: Not matrix based
    // this function converts the local model coordinates of the
    // sent render list into world coordinates, the results are stored
    // in the transformed vertex list (tvlist) within the renderlist

    // interate thru vertex list and transform all the model/local
    // coords to world coords by translating the vertex list by
    // the amount world_pos and storing the results in tvlist[]
    // is this polygon valid?

    if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            // all good, let's transform
            for (int vertex = 0; vertex < 3; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&curr_poly->vlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            } // end for vertex

        } // end for poly
    }     // end if local
    else  // TRANSFORM_TRANS_ONLY
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            // acquire current polygon
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            // transform this polygon if and only if it's not clipped, not culled,
            // active, and visible, note however the concept of "backface" is
            // irrelevant in a wire frame engine though
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue; // move onto next poly

            for (int vertex = 0; vertex < 3; vertex++)
            {
                // translate vertex
                VECTOR4D_Add(&curr_poly->tvlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            } // end for vertex

        } // end for poly

    } // end else

} // end Model_To_World_RENDERLIST4DV1

void Build_CAM4DV1_Matrix_Euler(CAM4DV1_PTR cam, int cam_rot_seq)
{
    // this creates a camera matrix based on Euler angles
    // and stores it in the sent camera object
    // if you recall from chapter 6 to create the camera matrix
    // we need to create a transformation matrix that looks like:

    // Mcam = mt(-1) * my(-1) * mx(-1) * mz(-1)
    // that is the inverse of the camera translation matrix mutilplied
    // by the inverses of yxz, in that order, however, the order of
    // the rotation matrices is really up to you, so we aren't going
    // to force any order, thus its programmable based on the value
    // of cam_rot_seq which can be any value CAM_ROT_SEQ_XYZ where
    // XYZ can be in any order, YXZ, ZXY, etc.

    MATRIX4X4 mt_inv, // inverse camera translation matrix
        mx_inv,       // inverse camera x axis rotation matrix
        my_inv,       // inverse camera y axis rotation matrix
        mz_inv,       // inverse camera z axis rotation matrix
        mrot,         // concatenated inverse rotation matrices
        mtmp;         // temporary working matrix

    // step 1: create the inverse translation matrix for the camera
    // position
    Mat_Init_4X4(&mt_inv, 1, 0, 0, 0,
                 0, 1, 0, 0,
                 0, 0, 1, 0,
                 -cam->pos.x, -cam->pos.y, -cam->pos.z, 1);

    // step 2: create the inverse rotation sequence for the camera
    // rember either the transpose of the normal rotation matrix or
    // plugging negative values into each of the rotations will result
    // in an inverse matrix

    // first compute all 3 rotation matrices

    // extract out euler angles
    float theta_x = cam->dir.x;
    float theta_y = cam->dir.y;
    float theta_z = cam->dir.z;

    // compute the sine and cosine of the angle x
    float cos_theta = Fast_Cos(theta_x);  // no change since cos(-x) = cos(x)
    float sin_theta = -Fast_Sin(theta_x); // sin(-x) = -sin(x)

    // set the matrix up
    Mat_Init_4X4(&mx_inv, 1, 0, 0, 0,
                 0, cos_theta, sin_theta, 0,
                 0, -sin_theta, cos_theta, 0,
                 0, 0, 0, 1);

    // compute the sine and cosine of the angle y
    cos_theta = Fast_Cos(theta_y);  // no change since cos(-x) = cos(x)
    sin_theta = -Fast_Sin(theta_y); // sin(-x) = -sin(x)

    // set the matrix up
    Mat_Init_4X4(&my_inv, cos_theta, 0, -sin_theta, 0,
                 0, 1, 0, 0,
                 sin_theta, 0, cos_theta, 0,
                 0, 0, 0, 1);

    // compute the sine and cosine of the angle z
    cos_theta = Fast_Cos(theta_z);  // no change since cos(-x) = cos(x)
    sin_theta = -Fast_Sin(theta_z); // sin(-x) = -sin(x)

    // set the matrix up
    Mat_Init_4X4(&mz_inv, cos_theta, sin_theta, 0, 0,
                 -sin_theta, cos_theta, 0, 0,
                 0, 0, 1, 0,
                 0, 0, 0, 1);

    // now compute inverse camera rotation sequence
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
    } // end switch

    // now mrot holds the concatenated product of inverse rotation matrices
    // multiply the inverse translation matrix against it and store in the
    // camera objects' camera transform matrix we are done!
    Mat_Mul_4X4(&mt_inv, &mrot, &cam->mcam);

} // end Build_CAM4DV1_Matrix_Euler

void World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
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
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue; // move onto next poly

        // all good, let's transform
        for (int vertex = 0; vertex < 3; vertex++)
        {
            // transform the vertex by the mcam matrix within the camera
            // it better be valid!
            POINT4D presult; // hold result of each transformation

            // transform point
            Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], &cam->mcam, &presult);

            // store result back
            VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
        } // end for vertex

    } // end for poly

} // end World_To_Camera_RENDERLIST4DV1

void Camera_To_Perspective_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
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
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
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

} // end Camera_To_Perspective_RENDERLIST4DV1

void Perspective_To_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam)
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
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        // is this polygon valid?
        // transform this polygon if and only if it's not clipped, not culled,
        // active, and visible, note however the concept of "backface" is
        // irrelevant in a wire frame engine though
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
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

} // end Perspective_To_Screen_RENDERLIST4DV1

void Mat_Mul_VECTOR4D_4X4(VECTOR4D_PTR va, MATRIX4X4_PTR mb, VECTOR4D_PTR vprod)
{
    // this function multiplies a VECTOR4D against a
    // 4x4 matrix - ma*mb and stores the result in mprod
    // the function makes no assumptions

    for (int col = 0; col < 4; col++)
    {
        // compute dot product from row of ma
        // and column of mb
        float sum = 0; // used to hold result

        for (int row = 0; row < 4; row++)
        {
            // add in next product pair
            sum += (va->M[row] * mb->M[row][col]);
        } // end for index

        // insert resulting col element
        vprod->M[col] = sum;

    } // end for col

} // end Mat_Mul_VECTOR4D_4X4

void VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb, VECTOR4D_PTR vsum)
{
    // this function adds va+vb and return it in vsum
    vsum->x = va->x + vb->x;
    vsum->y = va->y + vb->y;
    vsum->z = va->z + vb->z;
    vsum->w = 1;

} // end VECTOR4D_Add

////////////////////////////////////////////////////////////

VECTOR4D VECTOR4D_Add(VECTOR4D_PTR va, VECTOR4D_PTR vb)
{
    // this function adds va+vb and returns the result on
    // the stack
    VECTOR4D vsum;

    vsum.x = va->x + vb->x;
    vsum.y = va->y + vb->y;
    vsum.z = va->z + vb->z;
    vsum.w = 1;

    // return result
    return (vsum);

} // end VECTOR4D_Add