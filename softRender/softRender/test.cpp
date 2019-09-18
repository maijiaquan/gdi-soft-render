void Transform_RENDERLIST4DV1(RENDERLIST4DV1_PTR rend_list, MATRIX4X4_PTR mt, int coord_select)
{
    ...
    for (int poly = 0; poly < rend_list->num_polys; poly++)
    {
        POLYF4DV1_PTR curr_poly = rend_list->poly_ptrs[poly];
        for (int vertex = 0; vertex < 3; vertex++)
        {
            POINT4D point;
            Mat_Mul_VECTOR4D_4X4(&curr_poly->vlist[vertex], mt, &point);
            VECTOR4D_COPY(&curr_poly->vlist[vertex], &point);
        }
    }
    ...
}