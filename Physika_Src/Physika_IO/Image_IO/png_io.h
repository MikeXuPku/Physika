/*
 * @file png_io.h 
 * @Brief load/save png file
 * @author Fei Zhu
 * 
 * This file is part of Physika, a versatile physics simulation library.
 * Copyright (C) 2013 Physika Group.
 *
 * This Source Code Form is subject to the terms of the GNU General Public License v2.0. 
 * If a copy of the GPL was not distributed with this file, you can obtain one at:
 * http://www.gnu.org/licenses/gpl-2.0.html
 *
 */

#ifndef PHYSIKA_IO_IMAGE_IO_PNG_IO_H_
#define PHYSIKA_IO_IMAGE_IO_PNG_IO_H_

#include <string>

namespace Physika{

class PngIO
{
public:
    PngIO(){}
    ~PngIO(){}
    /* load PNG image from given file, return the image data in row order
     * if load fails, return NULL 
     * memory of the image data needs to be released by the caller
     */
    static unsigned char* load(const std::string &filename, int &width, int &height);

    /* save image data to file, the image data is in row order
     */
    static void save(const std::string &filename, int width, int height, const unsigned char *image_data);
protected:

};

} //end of namespace Physika

#endif //PHYSIKA_IO_IMAGE_IO_PNG_IO_H_










