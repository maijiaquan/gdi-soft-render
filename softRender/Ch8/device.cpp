#include "device.h"

// 画点
void device_pixel(device_t *device, int x, int y, IUINT32 color)
{
	if (((IUINT32)x) < (IUINT32)device->width && ((IUINT32)y) < (IUINT32)device->height)
	{
		device->framebuffer[y][x] = color;
	}
}


// 设备初始化，fb为外部帧缓存，非 NULL 将引用外部帧缓存（每行 4字节对齐）
void device_init(device_t *device, int width, int height, void *fb)
{
	int need = sizeof(void *) * (height * 2 + 1024) + width * height * 8;
	char *ptr = (char *)malloc(need + 64);
	char *framebuf, *zbuf;
	int j;
	assert(ptr);
	device->framebuffer = (IUINT32 **)ptr;
	device->zbuffer = (float **)(ptr + sizeof(void *) * height);
	ptr += sizeof(void *) * height * 2;
	device->texture = (IUINT32 **)ptr;
	ptr += sizeof(void *) * 1024;
	framebuf = (char *)ptr;
	zbuf = (char *)ptr + width * height * 4;
	ptr += width * height * 8;
	if (fb != NULL)
		framebuf = (char *)fb;
	for (j = 0; j < height; j++)
	{
		device->framebuffer[j] = (IUINT32 *)(framebuf + width * 4 * j);
		device->zbuffer[j] = (float *)(zbuf + width * 4 * j);
	}
	device->texture[0] = (IUINT32 *)ptr;
	device->texture[1] = (IUINT32 *)(ptr + 16);
	memset(device->texture[0], 0, 64);
	device->tex_width = 2;
	device->tex_height = 2;
	device->max_u = 1.0f;
	device->max_v = 1.0f;
	device->width = width;
	device->height = height;
	device->background = 0xc0c0c0;
	device->foreground = 0xc0c0c0;
	device->render_state = RENDER_STATE_WIREFRAME;
}

// 删除设备
void device_destroy(device_t *device)
{
	if (device->framebuffer)
		free(device->framebuffer);
	device->framebuffer = NULL;
	device->zbuffer = NULL;
	device->texture = NULL;
}

// 清空 framebuffer 和 zbuffer
void device_clear(device_t *device, int mode)
{
	int y, x, height = device->height;
	for (y = 0; y < device->height; y++)
	{
		IUINT32 *dst = device->framebuffer[y];
		// IUINT32 cc = (height - 1 - y) * 230 / (height - 1);
		// IUINT32 cc = (height - 1 - 0) * 255 / (height - 1);
		IUINT32 cc = (height - 1 - 0) * 0 / (height - 1); //黑色背景
		int R = 0;
		int G = 0;
		int B = 0;
		//cc = (cc << 16) | (cc << 8) | cc;
		cc = (R << 16) | (G << 8) | B;

		if (mode == 0)
			cc = device->background;
		for (x = device->width; x > 0; dst++, x--)
			dst[0] = cc;
	}
	for (y = 0; y < device->height; y++)
	{
		float *dst = device->zbuffer[y];
		for (x = device->width; x > 0; dst++, x--)
			dst[0] = 0.0f;
	}
}

// 绘制线段
void device_draw_line(device_t *device, int x1, int y1, int x2, int y2, IUINT32 c)
{
	// std::cout<<"x1 = "<<x1<<"y1 = "<<y1<<std::endl;
	// std::cout<<"x2 = "<<x2<<"y2 = "<<y2<<std::endl;

	int x, y, rem = 0;
	if (x1 == x2 && y1 == y2)
	{
		device_pixel(device, x1, y1, c);
	}
	else if (x1 == x2)
	{
		int inc = (y1 <= y2) ? 1 : -1;
		for (y = y1; y != y2; y += inc)
			device_pixel(device, x1, y, c);
		device_pixel(device, x2, y2, c);
	}
	else if (y1 == y2)
	{
		int inc = (x1 <= x2) ? 1 : -1;
		for (x = x1; x != x2; x += inc)
			device_pixel(device, x, y1, c);
		device_pixel(device, x2, y2, c);
	}
	else
	{
		int dx = (x1 < x2) ? x2 - x1 : x1 - x2;
		int dy = (y1 < y2) ? y2 - y1 : y1 - y2;
		if (dx >= dy)
		{
			if (x2 < x1)
				x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; x <= x2; x++)
			{
				device_pixel(device, x, y, c);
				rem += dy;
				if (rem >= dx)
				{
					rem -= dx;
					y += (y2 >= y1) ? 1 : -1;
					device_pixel(device, x, y, c);
				}
			}
			device_pixel(device, x2, y2, c);
		}
		else
		{
			if (y2 < y1)
				x = x1, y = y1, x1 = x2, y1 = y2, x2 = x, y2 = y;
			for (x = x1, y = y1; y <= y2; y++)
			{
				device_pixel(device, x, y, c);
				rem += dx;
				if (rem >= dy)
				{
					rem -= dy;
					x += (x2 >= x1) ? 1 : -1;
					device_pixel(device, x, y, c);
				}
			}
			device_pixel(device, x2, y2, c);
		}
	}
}

void DrawTopTriangle(device_t *device, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
{

    float dx_right, // the dx/dy ratio of the right edge of line
        dx_left,    // the dx/dy ratio of the left edge of line
        xs, xe,     // the starting and ending points of the edges
        height;     // the height of the triangle

    int temp_x, // used during sorting as temps
        temp_y,
        right, // used by clipping
        left;

    // test order of x1 and x2
    //保证 x1 < x2
    if (x2 < x1)
    {
        temp_x = x2;
        x2 = x1;
        x1 = temp_x;
    } // end if swap

    // compute delta's
    height = y3 - y1;

    dx_left = (x3 - x1) / height;
    dx_right = (x3 - x2) / height;

    // set starting points
    xs = (float)x1;
    xe = (float)x2; // +(float)0.5;

    // perform y clipping
    if (y1 < min_clip_y)
    {
        // compute new xs and ys
        xs = xs + dx_left * (float)(-y1 + min_clip_y);
        xe = xe + dx_right * (float)(-y1 + min_clip_y);

        // reset y1
        y1 = min_clip_y;

    } // end if top is off screen

    if (y3 > max_clip_y)
        y3 = max_clip_y;

    // compute starting address in video memory

    // test if x clipping is needed
    if (x1 >= min_clip_x && x1 <= max_clip_x &&
        x2 >= min_clip_x && x2 <= max_clip_x &&
        x3 >= min_clip_x && x3 <= max_clip_x)
    {
        // draw the triangle
        //for (temp_y = y1; temp_y <= y3; temp_y++, dest_addr += mempitch)
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // draw the line

            device_draw_line(device, xs, temp_y, xe, temp_y, color);

            xs += dx_left;
            xe += dx_right;

        } // end for

    } // end if no x clipping needed
    else
    {
        // clip x axis with slower version

        // draw the triangle
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // do x clip
            left = (int)xs;
            right = (int)xe;

            // adjust starting point and ending point
            xs += dx_left;
            xe += dx_right;

            // clip line
            if (left < min_clip_x)
            {
                left = min_clip_x;

                if (right < min_clip_x)
                    continue;
            }

            if (right > max_clip_x)
            {
                right = max_clip_x;

                if (left > max_clip_x)
                    continue;
            }

            // draw the line
            // IUINT32 c = (0 << 16) | (255 << 8) | 0;
            device_draw_line(device, left, temp_y, right, temp_y, color);
        } // end for

    } // end else x clipping needed
}

void DrawDownTriangle(device_t *device, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
{

    float dx_right, // the dx/dy ratio of the right edge of line
        dx_left,    // the dx/dy ratio of the left edge of line
        xs, xe,     // the starting and ending points of the edges
        height;     // the height of the triangle

    int temp_x, // used during sorting as temps
        temp_y,
        right, // used by clipping
        left;

    // test order of x1 and x2
    if (x3 < x2)
    {
        temp_x = x2;
        x2 = x3;
        x3 = temp_x;
    } // end if swap

    // compute delta's
    height = y3 - y1;

    dx_left = (x2 - x1) / height;
    dx_right = (x3 - x1) / height;

    // set starting points
    xs = (float)x1;
    xe = (float)x1; // +(float)0.5;

    // perform y clipping
    if (y1 < min_clip_y)
    {
        // compute new xs and ys
        xs = xs + dx_left * (float)(-y1 + min_clip_y);
        xe = xe + dx_right * (float)(-y1 + min_clip_y);

        // reset y1
        y1 = min_clip_y;

    } // end if top is off screen

    if (y3 > max_clip_y)
        y3 = max_clip_y;

    // compute starting address in video memory

    // test if x clipping is needed
    if (x1 >= min_clip_x && x1 <= max_clip_x &&
        x2 >= min_clip_x && x2 <= max_clip_x &&
        x3 >= min_clip_x && x3 <= max_clip_x)
    {
        // draw the triangle
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // draw the line
            device_draw_line(device, xs, temp_y, xe, temp_y, color);

            // adjust starting point and ending point
            xs += dx_left;
            xe += dx_right;

        } // end for

    } // end if no x clipping needed
    else
    {
        // clip x axis with slower version

        // draw the triangle
        for (temp_y = y1; temp_y <= y3; temp_y++)
        {
            // do x clip
            left = (int)xs;
            right = (int)xe;

            // adjust starting point and ending point
            xs += dx_left;
            xe += dx_right;

            // clip line
            if (left < min_clip_x)
            {
                left = min_clip_x;

                if (right < min_clip_x)
                    continue;
            }

            if (right > max_clip_x)
            {
                right = max_clip_x;

                if (left > max_clip_x)
                    continue;
            }
            // draw the line
            device_draw_line(device, left, temp_y, right, temp_y, color);
        } // end for

    } // end else x clipping needed
}

void DrawTrianglePureColor(device_t *device, int x1, int y1, int x2, int y2, int x3, int y3, IUINT32 color)
{

		int temp_x, // used for sorting
			temp_y,
			new_x;

		// test for h lines and v lines
		if ((x1 == x2 && x2 == x3) || (y1 == y2 && y2 == y3))
			return;

		//根据y的大小排序
		if (y2 < y1)
		{
			temp_x = x2;
			temp_y = y2;
			x2 = x1;
			y2 = y1;
			x1 = temp_x;
			y1 = temp_y;
		} // end if

		// now we know that p1 and p2 are in order
		if (y3 < y1)
		{
			temp_x = x3;
			temp_y = y3;
			x3 = x1;
			y3 = y1;
			x1 = temp_x;
			y1 = temp_y;
		} // end if

		// finally test y3 against y2
		if (y3 < y2)
		{
			temp_x = x3;
			temp_y = y3;
			x3 = x2;
			y3 = y2;
			x2 = temp_x;
			y2 = temp_y;

		} // end if

		// do trivial rejection tests for clipping
		if (y3 < min_clip_y || y1 > max_clip_y ||
			(x1 < min_clip_x && x2 < min_clip_x && x3 < min_clip_x) ||
			(x1 > max_clip_x && x2 > max_clip_x && x3 > max_clip_x))
			return;

		// test if top of triangle is flat
		if (y1 == y2)
		{

			DrawTopTriangle(device, x1, y1, x2, y2, x3, y3, color);
		} // end if
		else if (y2 == y3)
		{
			DrawDownTriangle(device, x1, y1, x2, y2, x3, y3, color);
		} // end if bottom is flat
		else
		{
			// draw each sub-triangle
			new_x = x1 + (int)(0.5 + (float)(y2 - y1) * (float)(x3 - x1) / (float)(y3 - y1));
			DrawDownTriangle(device, x1, y1, new_x, y2, x2, y2, color);
			DrawTopTriangle(device, x2, y2, new_x, y2, x3, y3, color);
		} // end else

}