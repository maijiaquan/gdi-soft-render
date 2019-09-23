#include "transform.h"

void Reset_RENDERLIST4DV1(RENDERLIST4DV1_PTR renderList)
{
    renderList->num_polys = 0;
}

void Reset_OBJECT4DV1(OBJECT4DV1_PTR obj)
{
    RESET_BIT(obj->state, OBJECT4DV1_STATE_CULLED);

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        POLY4DV1_PTR curr_poly = &obj->plist[poly];

        if (!(curr_poly->state & POLY4DV1_STATE_ACTIVE))
            continue;
        RESET_BIT(curr_poly->state, POLY4DV1_STATE_CLIPPED);
        RESET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
    }
}

//����Ⱦ�б��ÿһ������ν��о���任
void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, MATRIX4X4_PTR mt, int coord_select)
{
    switch (coord_select)
    {
    case TRANSFORM_LOCAL_ONLY:
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue;

            for (int vertex = 0; vertex < 3; vertex++)
            {
                POINT4D point;
                Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &point);
                VECTOR4D_COPY(&curr_poly->vlist[vertex], &point);
            }
        }
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

//ģ�����굽���������ת��
void Model_To_World_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, POINT4D_PTR world_pos, int coord_select)
{
    if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
                (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
                (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue;

            for (int vertex = 0; vertex < 3; vertex++)
            {
                VECTOR4D_Add(&curr_poly->vlist[vertex], world_pos, &curr_poly->tvlist[vertex]); //ֱ�ӽ��������ӷ�
            }
        }
    }
    else
    {
        for (int poly = 0; poly < rend_list->num_polys; poly++)
        {
            POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
            if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) || (curr_poly->state & POLY4DV1_STATE_CLIPPED) || (curr_poly->state & POLY4DV1_STATE_BACKFACE))
                continue;

            for (int vertex = 0; vertex < 3; vertex++)
            {
                VECTOR4D_Add(&curr_poly->tvlist[vertex], world_pos, &curr_poly->tvlist[vertex]);
            }
        }
    }
}

//�������굽��������ת��
void World_To_Camera_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list,
                                    CAM4DV1_PTR cam)
{
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];

        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) || (curr_poly->state & POLY4DV1_STATE_CLIPPED) || (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue;

        for (int vertex = 0; vertex < 3; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&curr_poly->tvlist[vertex], &cam->mcam, &presult);
            VECTOR4D_COPY(&curr_poly->tvlist[vertex], &presult);
        }
    }
}

//������굽͸������
void Camera_To_Perspective_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam)
{
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) || (curr_poly->state & POLY4DV1_STATE_CLIPPED) || (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue;

        for (int vertex = 0; vertex < 3; vertex++)
        {
            float z = curr_poly->tvlist[vertex].z;
            curr_poly->tvlist[vertex].x = cam->view_dist * curr_poly->tvlist[vertex].x / z;
            curr_poly->tvlist[vertex].y = cam->view_dist * curr_poly->tvlist[vertex].y * cam->aspect_ratio / z;
        }
    }
}

//͸�����굽��Ļ����
void Perspective_To_Screen_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, CAM4DV1_PTR cam)
{
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        if ((curr_poly == NULL) || !(curr_poly->state & POLY4DV1_STATE_ACTIVE) ||
            (curr_poly->state & POLY4DV1_STATE_CLIPPED) ||
            (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue;

        float alpha = (0.5 * cam->viewport_width - 0.5);
        float beta = (0.5 * cam->viewport_height - 0.5);

        for (int vertex = 0; vertex < 3; vertex++)
        {
            curr_poly->tvlist[vertex].x = alpha + alpha * curr_poly->tvlist[vertex].x;
            curr_poly->tvlist[vertex].y = beta - beta * curr_poly->tvlist[vertex].y;
        }
    }
}

//////////////////////////////// ���� //////////////////////////////////////////

//ͨ�ñ任����
void Transform_OBJECT4DV1(OBJECT4DV1_PTR obj, MATRIX4X4_PTR mt, int coord_select, int transform_basis)
{
    switch (coord_select)
    {
    case TRANSFORM_LOCAL_ONLY: //local��local
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex], mt, &presult);
            VECTOR4D_COPY(&obj->vlist_local[vertex], &presult);
        }
    }
    break;

    case TRANSFORM_TRANS_ONLY: //trans��trans
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex], mt, &presult);
            VECTOR4D_COPY(&obj->vlist_trans[vertex], &presult);
        }
    }
    break;

    case TRANSFORM_LOCAL_TO_TRANS: //local��trans
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            POINT4D presult;
            Mat_Mul_VECTOR4D_4X4(&obj->vlist_local[vertex], mt, &obj->vlist_trans[vertex]);
        }
    }
    break;

    default:
        break;
    }

    if (transform_basis)
    {
        VECTOR4D vresult;

        // rotate ux of basis
        Mat_Mul_VECTOR4D_4X4(&obj->ux, mt, &vresult);
        VECTOR4D_COPY(&obj->ux, &vresult);

        // rotate uy of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uy, mt, &vresult);
        VECTOR4D_COPY(&obj->uy, &vresult);

        // rotate uz of basis
        Mat_Mul_VECTOR4D_4X4(&obj->uz, mt, &vresult);
        VECTOR4D_COPY(&obj->uz, &vresult);
    }
}

//��������ı任�����������obj->vlist_trans
void Model_To_World_OBJECT4DV1(OBJECT4DV1_PTR obj, int coord_select)
{
    if (coord_select == TRANSFORM_LOCAL_TO_TRANS)
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            VECTOR4D_Add(&obj->vlist_local[vertex], &obj->world_pos, &obj->vlist_trans[vertex]);
        }
    }
    else
    {
        for (int vertex = 0; vertex < obj->num_vertices; vertex++)
        {
            VECTOR4D_Add(&obj->vlist_trans[vertex], &obj->world_pos, &obj->vlist_trans[vertex]);
        }
    }
}

//������obj->vlist_trans��ִ������任�����������obj->vlist_trans
void World_To_Camera_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        POINT4D presult;
        Mat_Mul_VECTOR4D_4X4(&obj->vlist_trans[vertex], &cam->mcam, &presult);
        VECTOR4D_COPY(&obj->vlist_trans[vertex], &presult);
    }
}

//������obj->vlist_trans��ִ��͸�ӱ任�����������obj->vlist_trans
void Camera_To_Perspective_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{
    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        float z = obj->vlist_trans[vertex].z;
        obj->vlist_trans[vertex].x = cam->view_dist * obj->vlist_trans[vertex].x / z;
        obj->vlist_trans[vertex].y = cam->view_dist * obj->vlist_trans[vertex].y * cam->aspect_ratio / z;
    }
}

//������obj->vlist_trans��ִ����Ļ�任�����������obj->vlist_trans
void Perspective_To_Screen_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{

    float alpha = (0.5 * cam->viewport_width - 0.5);
    float beta = (0.5 * cam->viewport_height - 0.5);

    for (int vertex = 0; vertex < obj->num_vertices; vertex++)
    {
        obj->vlist_trans[vertex].x = alpha + alpha * obj->vlist_trans[vertex].x;
        obj->vlist_trans[vertex].y = beta - beta * obj->vlist_trans[vertex].y;
    }
}


void Remove_Backfaces_OBJECT4DV1(OBJECT4DV1_PTR obj, CAM4DV1_PTR cam)
{
    if (obj->state & OBJECT4DV1_STATE_CULLED)
        return;

    for (int poly = 0; poly < obj->num_polys; poly++)
    {
        POLY4DV1_PTR curr_poly = &obj->plist[poly];

        if (!(curr_poly->state & POLY4DV1_STATE_ACTIVE) || (curr_poly->state & POLY4DV1_STATE_CLIPPED) || (curr_poly->attr & POLY4DV1_ATTR_2SIDED) || (curr_poly->state & POLY4DV1_STATE_BACKFACE))
            continue; 

        int vindex_0 = curr_poly->vert[0];
        int vindex_1 = curr_poly->vert[1];
        int vindex_2 = curr_poly->vert[2];


        VECTOR4D u, v, n;

        // ��������u��v
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &obj->vlist_trans[vindex_1], &u);
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &obj->vlist_trans[vindex_2], &v);

        // ���㷨���� 
        VECTOR4D_Cross(&u, &v, &n);

        //����۲�����
        VECTOR4D view;
        VECTOR4D_Build(&obj->vlist_trans[vindex_0], &cam->pos, &view);

        float dp = VECTOR4D_Dot(&n, &view);

        if (dp <= 0.0)
            SET_BIT(curr_poly->state, POLY4DV1_STATE_BACKFACE);
    } 
}