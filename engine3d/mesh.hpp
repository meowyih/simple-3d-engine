#ifndef OCTILLION_ENGINE3D_MESH_HEADER
#define OCTILLION_ENGINE3D_MESH_HEADER

#include <iostream>
#include <vector>
#include <string>

#include "../jsonw/jsonw.hpp"
#include "face.hpp"

class Mesh
{
public:
  Mesh() {}
  void init(std::string jsonfile)
  {
    // reset data
    valid_ = false;
    faces_.clear();

    // open file
    std::ifstream fin(jsonfile);
    if (!fin.good())
      return;

    // read file into json
    octillion::JsonW json(fin);
    if (json.valid() == false)
      return;

    // get vertices array
    octillion::JsonW* vertices = json.get("vertices");
    if (vertices->size() == 0 || vertices->size() % 9 != 0)
      return;

    // feed data into faces_
    faces_.reserve(vertices->size() / 3);
    for (size_t idx = 0; idx < vertices->size(); idx = idx + 9)
    {
      faces_.push_back(
        Face(
          PointF( 
            vertices->get(idx)->frac(), 
            vertices->get(idx + 1)->frac(), 
            vertices->get(idx + 2)->frac() ),
          PointF(
            vertices->get(idx + 3)->frac(),
            vertices->get(idx + 4)->frac(),
            vertices->get(idx + 5)->frac()),
          PointF(
            vertices->get(idx + 6)->frac(),
            vertices->get(idx + 7)->frac(),
            vertices->get(idx + 8)->frac())
        )
      );
    }

    valid_ = true;
  }

public:
  bool valid_ = false;
  std::vector<Face> faces_;
};


#endif
