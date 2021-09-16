// Copyright (c) 2017  GeometryFactory Sarl (France)
// All rights reserved.
//
// This file is part of CGAL (www.cgal.org).
// You can redistribute it and/or modify it under the terms of the GNU
// General Public License as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// Licensees holding a valid commercial license may use this file in
// accordance with the commercial license agreement provided with the software.
//
// This file is provided AS IS with NO WARRANTY OF ANY KIND, INCLUDING THE
// WARRANTY OF DESIGN, MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE.
//
// $URL: https://github.com/CGAL/cgal/blob/releases/CGAL-4.14.3/Three/include/CGAL/Three/Primitive_container.h $
// $Id: Primitive_container.h 5b4cf63 2018-06-26T14:57:44+02:00 Maxime Gimeno
// SPDX-License-Identifier: GPL-3.0+
//
// Author(s)     : Maxime Gimeno

#ifndef PRIMITIVE_CONTAINER_H
#define PRIMITIVE_CONTAINER_H

#include <CGAL/license/Three.h>


#include <CGAL/Three/Buffer_objects.h>
#include <CGAL/Three/Viewer_interface.h>
#include <CGAL/Three/Scene_item_rendering_helper.h>

using namespace CGAL::Three;

#ifdef demo_framework_EXPORTS
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_EXPORT
#else
#  define DEMO_FRAMEWORK_EXPORT Q_DECL_IMPORT
#endif
class QOpenGLFramebufferObject;
struct D;
namespace CGAL {
namespace Three {


//!
//! \brief The Primitive_container struct provides a base for the OpenGL data wrappers.
//!
class DEMO_FRAMEWORK_EXPORT Primitive_container
{
public:

  //!
  //! \brief Primitive_container ructor.
  //! \param program the `QOpenGLShaderProgram` used by the VAOs.
  //! \param indexed must be `true` if the data is indexed, `false` otherwise.
  //!
  Primitive_container(int program, bool indexed);

  virtual ~Primitive_container();

  /*!
     * \brief bindUniformValues sets the uniform variables for the concerned shaders.
     *
     * Such variables are valid at every step of the pipeline. For example,
     * the ModelViewProjection matrix, the uniform color or the is_selected state are uniform values.
     * This function is called in the `draw()`function.
     * \attention `Vbo`s data should be allocated for this function to be effective.
     * \attention This should only be called once the `Vao`s and the `Vbo`s are created, in a
     * valid OpenGL context.
     * \param viewer the active `Viewer_interface`
     */
  void bindUniformValues(CGAL::Three::Viewer_interface* viewer) ;

  //!
  //! \brief draw is the function that actually renders the data.
  //! \param viewer the active `Viewer_interface`.
  //! \param is_color_uniform should be `true` if the item is unicolor.
  //!
  virtual void draw(CGAL::Three::Viewer_interface* viewer,
                    bool is_color_uniform)  = 0;

  //!
  //! \brief initializeBuffers sends the data to the GPU memory.
  //!
  //! It actually fills up the buffers with the data provided by `Vbo::allocate()`;
  //! \param viewer the active `Viewer_interface`.
  //!
  virtual void initializeBuffers(CGAL::Three::Viewer_interface* viewer) ;

  //!
  //! \brief initGL initializes the OpenGL containers.
  //! \attention It must be called within a valid OpenGL context. The `draw()` function of an item is always a safe place to call this.
  //!
  //! \param viewer the active `Viewer_interface`.
  virtual void initGL(CGAL::Three::Viewer_interface* viewer)  = 0;

  //!
  //! \brief removeViewer deletes and removes the Vao assigned to `viewer` from `Vaos`.
  //! \param viewer the `Viewer_interface` to remove.
  //!
  void removeViewer(CGAL::Three::Viewer_interface* viewer) ;

  //!
  //! \brief reset_vbos de-allocates the `Vbo`s. It must be called when the `Vbo`s data is updated.
  //!
  void reset_vbos(Scene_item_rendering_helper::Gl_data_names);

  //!\todo is it better to have a non  draw() or to have 99% of the setters  ?

  //!
  //! \brief setFlatDataSize sets the number of un-indexed
  //! vertices of this container.
  //!
  //! If the vertices are indexed, you can ignore this function.
  //!
  void setFlatDataSize(std::size_t);
  //!
  //! \brief setIdxDataSize sets the number of indexed
  //! vertices of this container.
  //!
  //! If the vertices are not indexed, you can ignore this function.
  //!
  void setIdxSize(std::size_t);
  //!
  //! \brief setCenterSize sets the number of instances of
  //! the item in this container.
  //!
  //! If the program of this container is not instanced, you can ignore this function.
  //!
  void setCenterSize(std::size_t);

  //! Returns `true` if the container `Vao`s and `Vbo`s are
  //! created in the context of `viewer`.
  bool isGLInit(Viewer_interface* viewer) const;
  //!
  //! \brief allocate sets the data for a `Vbo`.
  //! \param vbo_id the index of the `Vbo` in this container vector.
  //! \param data the data to give to the `Vbo`.
  //! \param datasize the size in bytes of `data`.
  //!
  void allocate(std::size_t vbo_id, void* data, int datasize);
  //!
  //! \name Setters for the shaders parameters.
  //!@{
  
  //! Setter for the "selected" uniform parameter.
  void setSelected(bool);
  //! Setter for the "color" parameter.
  void setColor(QColor);
  //!@}
protected:
  //!
  //! \brief setVao sets the `Vao` corresponding to `viewer` of this container.
  //!
  void setVao(Viewer_interface* viewer, Vao*) ;
  //!
  //! \brief setVbos sets the vector of `Vbo`s for this container.
  //! It must contain at least all the `Vbo`s that will be used by the `Vao`s.
  //!
  void setVbos(std::vector<Vbo*>);
  //!
  //! \brief setVbo sets the `vbo_id`th `Vbo` of this container to `vbo`.
  //! \param vbo_id
  //! \param vbo
  //!
  void setVbo(std::size_t vbo_id, Vbo* vbo);
  //!
  //! Use this to specify if the container `Vao`s and `Vbo`s are
  //! created in the context of `viewer`.
  //!
  void setGLInit(Viewer_interface* viewer, bool) ;
  //!
  //! \brief Returns the `Vao` bound to `viewer`.
  //!
  Vao* getVao(Viewer_interface* viewer)const;
  //!
  //! \brief getVbo returns the `id`th Vbo of this container.
  //!
  Vbo *getVbo(std::size_t id)const;
  //!
  //! \brief getProgram returns the `OpenGL_program_IDs` used with this container.
  //!
  int getProgram()const;
  //!
  //! \brief isDataIndexed specifies if the data is indexed or not. This matters for the internal drawing functions.
  //!
  bool isDataIndexed();

  //!
  //!Use this to specify if the container `Vbo`s are filled for `viewer`.
  //!
  void setInit(Viewer_interface* viewer, bool) ;
  //!
  //! \brief isInit returns `true` if the container `Vbo`s are filled for `viewer`.
  bool isInit(Viewer_interface* viewer) const;
  //! \brief getFlatDataSize returns the number of un-indexed
  //! vertices.
  std::size_t getFlatDataSize()const;
  //! \brief getIdxSize returns the number of indexed
  //! vertices.
  std::size_t getIdxSize()const;
  //! \brief getCenterSize returns the number of instances of
  //! the item in this container.
  std::size_t getCenterSize()const;
  //! \name Getters for the shaders parameters.
  //!@{
  
  //! getter for the "selected" parameter
  bool isSelected()const;
  //! getter for the "color" parameter
  QColor getColor()const;
  //! @}
  //!
private:
  friend struct D;
  mutable D* d;
}; //end of class Triangle_container

}
}

#endif // PRIMITIVE_CONTAINER_H
