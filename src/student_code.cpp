#include "student_code.h"
#include "mutablePriorityQueue.h"

using namespace std;

namespace CGL
{

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (class member).
   *
   * @param points A vector of points in 2D
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector2D> BezierCurve::evaluateStep(std::vector<Vector2D> const &points)
  { 
    std::vector<Vector2D> results;
    // For each set of 2 points:
    for (int i = 0; i < points.size() - 1; i++) {
      // Create a new point and set its value to the interp of the two vectors
      float resultX = points[i].x * (1 - t) + points[i + 1].x * t;
      float resultY = points[i].y * (1 - t) + points[i + 1].y * t;
      
      Vector2D thisStep(resultX, resultY);
      
      results.push_back(thisStep);
    }
    // Stop at 1 less than the number of vectors
    
    return results;
  }

  /**
   * Evaluates one step of the de Casteljau's algorithm using the given points and
   * the scalar parameter t (function parameter).
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return A vector containing intermediate points or the final interpolated vector
   */
  std::vector<Vector3D> BezierPatch::evaluateStep(std::vector<Vector3D> const &points, double t) const
  {
    std::vector<Vector3D> results;
    // For each set of 2 points:
    for (int i = 0; i < points.size() - 1; i++) {
      // Create a new point and set its value to the interp of the two vectors
      float resultX = points[i].x * (1 - t) + points[i + 1].x * t;
      float resultY = points[i].y * (1 - t) + points[i + 1].y * t;
      float resultZ = points[i].z * (1 - t) + points[i + 1].z * t;

      Vector3D thisStep(resultX, resultY, resultZ);

      results.push_back(thisStep);
    }
    // Stop at 1 less than the number of vectors

    return results;
  }

  /**
   * Fully evaluates de Casteljau's algorithm for a vector of points at scalar parameter t
   *
   * @param points    A vector of points in 3D
   * @param t         Scalar interpolation parameter
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate1D(std::vector<Vector3D> const &points, double t) const
  {
    vector<Vector3D> intermed = points;
    while (intermed.size() > 1) {
      intermed = evaluateStep(intermed, t);
    }
    return intermed.at(0);
  }

  /**
   * Evaluates the Bezier patch at parameter (u, v)
   *
   * @param u         Scalar interpolation parameter
   * @param v         Scalar interpolation parameter (along the other axis)
   * @return Final interpolated vector
   */
  Vector3D BezierPatch::evaluate(double u, double v) const 
  {  
    vector<Vector3D> column;
    for (int i = 0; i < controlPoints.size(); i++) {
      column.push_back(evaluate1D(controlPoints.at(i), u));
    }

    return evaluate1D(column, v);
  }

  Vector3D Vertex::normal( void ) const
  {
    // TODO Part 3.
    // Returns an approximate unit normal at this vertex, computed by
    // taking the area-weighted average of the normals of neighboring
    // triangles, then normalizing.

    // Find each face's area by getting two vertices at a time.
    HalfedgeCIter h = halfedge();
    Vector3D normalTotal(0,0,0);
    do {
      VertexCIter neigh_0 = h->vertex();

      // Get the first neigboring vertex.
      HalfedgeCIter h_twin = h->twin();
      VertexCIter neigh_1 = h_twin->vertex();

      // Get the second neighboring vertex.
      HalfedgeCIter next_h = h_twin->next();
      VertexCIter neigh_2 = next_h->twin()->vertex();

      // Calculate the area using three vertices.
      Vector3D normal = cross(neigh_1->position, neigh_2->position) * -1;
      normalTotal = normalTotal + normal;

      h = h_twin->next();

    } while (h != halfedge());

    normalTotal.normalize();
    return normalTotal;
  }

  EdgeIter HalfedgeMesh::flipEdge( EdgeIter e0 )
  {
    // TODO Part 4.
    // This method should flip the given edge and return an iterator to the flipped edge.

    // Using the provided guide at http://15462.courses.cs.cmu.edu/fall2015content/misc/HalfedgeEdgeOpImplementationGuide.pdf

    // ## Phase 1: We collect all the relevant elements
    // Halfedges:
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h5->twin();

    // Vertices
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h2->vertex();
    VertexIter v3 = h5->vertex();

    // Edges
    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    // Faces
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    // ## Phase 1.5: We check if we are on a boundary
    if (f0->isBoundary() || f1->isBoundary()) {
      return e0;
    }

    // ## Phase 2: We reassign elements
    // Halfedges
    h0->setNeighbors(h1 , h3 , v3 , e0 , f0 );
    h1->setNeighbors(h2 , h7 , v2 , e2 , f0 );
    h2->setNeighbors(h0 , h8 , v0 , e3 , f0 );
    h3->setNeighbors(h4 , h0 , v2 , e0 , f1 );
    h4->setNeighbors(h5 , h9 , v3 , e4 , f1 );
    h5->setNeighbors(h3 , h6 , v1 , e1 , f1 );
    h6->setNeighbors(h6->next(), h5, v2, e1, h6->face());
    h7->setNeighbors(h7->next(), h1, v0, e2, h7->face());
    h8->setNeighbors(h8->next(), h2, v3, e3, h8->face());
    h9->setNeighbors(h9->next(), h4, v1, e4, h9->face());

    // Vertices
    v0->halfedge() = h2;
    v1->halfedge() = h5;
    v2->halfedge() = h3;
    v3->halfedge() = h0;

    // Edges
    e0->halfedge() = h0;
    e1->halfedge() = h5;
    e2->halfedge() = h1;
    e3->halfedge() = h2;
    e4->halfedge() = h4;

    // Faces
    f0->halfedge() = h0;
    f1->halfedge() = h3;

    // Done!
    return e0;
  }

  VertexIter HalfedgeMesh::splitEdge( EdgeIter e0 )
  {
    // TODO Part 5.
    // This method should split the given edge and return an iterator to the newly inserted vertex.
    // The halfedge of this vertex should point along the edge that was split, rather than the new edges.

    // ## Phase 1: Record existing elements
    // Halfedges:
    HalfedgeIter h0 = e0->halfedge();
    HalfedgeIter h1 = h0->next();
    HalfedgeIter h2 = h1->next();
    HalfedgeIter h3 = h0->twin();
    HalfedgeIter h4 = h3->next();
    HalfedgeIter h5 = h4->next();
    HalfedgeIter h6 = h1->twin();
    HalfedgeIter h7 = h2->twin();
    HalfedgeIter h8 = h4->twin();
    HalfedgeIter h9 = h5->twin();

    // Vertices
    VertexIter v0 = h0->vertex();
    VertexIter v1 = h3->vertex();
    VertexIter v2 = h2->vertex();
    VertexIter v3 = h5->vertex();

    // Edges
    EdgeIter e1 = h1->edge();
    EdgeIter e2 = h2->edge();
    EdgeIter e3 = h4->edge();
    EdgeIter e4 = h5->edge();

    // Faces
    FaceIter f0 = h0->face();
    FaceIter f1 = h3->face();

    // ## Phase 1.5: Stop if on a boundary
    if (f0->isBoundary() || f1->isBoundary()) {
      return VertexIter();
    }

    // ## Phase 1.75: Create new elements
    VertexIter v4 = newVertex();
    // Let's set the vertex's position.
    v4->position = (v0->position + v1->position) / 2;

    EdgeIter e5 = newEdge();
    EdgeIter e6 = newEdge();
    EdgeIter e7 = newEdge();
    HalfedgeIter h10 = newHalfedge();
    HalfedgeIter h11 = newHalfedge();
    HalfedgeIter h12 = newHalfedge();
    HalfedgeIter h13 = newHalfedge();
    HalfedgeIter h14 = newHalfedge();
    HalfedgeIter h15 = newHalfedge();
    FaceIter f2 = newFace();
    FaceIter f3 = newFace();

    // ## Phase 2: Set EVERYTHING AAAAAAAAAA SO MANY THINGs
    // Halfedges
    h0->setNeighbors(h1 , h3 , v0 , e0 , f0 );
    h1->setNeighbors(h2 , h10 , v4 , e5 , f0 );
    h2->setNeighbors(h0 , h7 , v2 , e2 , f0 );
    h3->setNeighbors(h4 , h0 , v4 , e0 , f1 );
    h4->setNeighbors(h5 , h8 , v0 , e3 , f1 );
    h5->setNeighbors(h3 , h13 , v3 , e7 , f1 );
    h6->setNeighbors(h6->next(), h12 , v2 , e1 , h6->face() );
    h7->setNeighbors(h7->next(), h2 , v0 , e2 , h7->face() );
    h8->setNeighbors(h8->next(), h4 , v3 , e3 , h8->face() );
    h9->setNeighbors(h9->next(), h14 , v1 , e4 , h9->face() );
    h10->setNeighbors(h11 , h1 , v2 , e5 , f2 );
    h11->setNeighbors(h12 , h15 , v4 , e6 , f2 );
    h12->setNeighbors(h10 , h6 , v1 , e1 , f2 );
    h13->setNeighbors(h14 , h5 , v4 , e7 , f3 );
    h14->setNeighbors(h15 , h9 , v3 , e4 , f3 );
    h15->setNeighbors(h13 , h11 , v1 , e6 , f3 );

    // Vertices
    v0->halfedge() = h0;
    v1->halfedge() = h15;
    v2->halfedge() = h10;
    v3->halfedge() = h5;
    v4->halfedge() = h3;

    // Edges
    e0->halfedge() = h0;
    e1->halfedge() = h12;
    e2->halfedge() = h2;
    e3->halfedge() = h4;
    e4->halfedge() = h14;
    e5->halfedge() = h1;
    e6->halfedge() = h11;
    e7->halfedge() = h5;

    // Faces
    f0->halfedge() = h0;
    f1->halfedge() = h3;
    f2->halfedge() = h10;
    f3->halfedge() = h13;

    return v4;
  }



  void MeshResampler::upsample( HalfedgeMesh& mesh )
  {
    // TODO Part 6.
    // This routine should increase the number of triangles in the mesh using Loop subdivision.
    // One possible solution is to break up the method as listed below.

    // 1. Compute new positions for all the vertices in the input mesh, using the Loop subdivision rule,
    // and store them in Vertex::newPosition. At this point, we also want to mark each vertex as being
    // a vertex of the original mesh.
    VertexIter activeVertex = mesh.verticesBegin();

    for (VertexIter activeVertex = mesh.verticesBegin(); activeVertex != mesh.verticesEnd(); activeVertex++) {
      // This is an old vertex.
      activeVertex->isNew = false;
      
      // This is an old vertex, so use the old vertex rule. Get n and u.
      float n = activeVertex->degree();
      float u = 3.0/16.0; // default for n=3
      if (n != 3) {
        u = 3.0 / (8.0 * n);
      }

      // Add up neighbor positions.
      Vector3D original_neighbor_position_sum(0, 0, 0);
      HalfedgeCIter h = activeVertex->halfedge();
      do {
        HalfedgeCIter h_twin = h->twin();
        VertexCIter v = h_twin->vertex();
        original_neighbor_position_sum = original_neighbor_position_sum + v->position;
        h = h_twin->next();
      } while (h != activeVertex->halfedge());

      // Set the newPosition.
      activeVertex->newPosition = (1 - n * u) * activeVertex->position + u * original_neighbor_position_sum;

    } while (activeVertex != mesh.verticesBegin());

    // 2. Compute the updated vertex positions associated with edges, and store it in Edge::newPosition.
    for (EdgeIter activeEdge = mesh.edgesBegin(); activeEdge != mesh.edgesEnd(); activeEdge++) {
      // Remember that it's an old edge.
      activeEdge->isNew = false;

      // The midpoint will become a new vertex, so use the new vertex equation.
      // Get A, B, C, D, where AB are on the edge and CD are on the touching faces.
      VertexIter a = activeEdge->halfedge()->vertex();
      VertexIter b = activeEdge->halfedge()->twin()->vertex();
      VertexIter c = activeEdge->halfedge()->next()->next()->vertex();
      VertexIter d = activeEdge->halfedge()->twin()->next()->next()->vertex();

      // Use the formula to set the new position.
      activeEdge->newPosition = (3.0 / 8.0) * (a->position + b->position) + (1.0 / 8.0) * (c->position + d->position);
    }
    
    // 3. Split every edge in the mesh, in any order. For future reference, we're also going to store some
    // information about which subdivide edges come from splitting an edge in the original mesh, and which edges
    // are new, by setting the flat Edge::isNew. Note that in this loop, we only want to iterate over edges of
    // the original mesh---otherwise, we'll end up splitting edges that we just split (and the loop will never end!)
    EdgeIter e = mesh.edgesBegin();
    while (e != mesh.edgesEnd()) {
      EdgeIter nextEdge = e;
      nextEdge++;

      // Let's check if this edge should be split. If either vertex is new, don't split it.
      VertexIter v0 = e->halfedge()->vertex();
      VertexIter v1 = e->halfedge()->twin()->vertex();
      bool shouldSplit = !v0->isNew && !v1->isNew;

      if (shouldSplit) {
        VertexIter newVertex = mesh.splitEdge(e);
        newVertex->isNew = true;

        // Record the new and old edges. Treat the split old edge as still old, even though it's technically two edges now.
        newVertex->halfedge()->edge()->isNew = false;
        newVertex->halfedge()->next()->next()->edge()->isNew = true;
        newVertex->halfedge()->next()->next()->twin()->next()->next()->edge()->isNew = false;
        newVertex->halfedge()->twin()->next()->edge()->isNew = true;

        // Copy over the new position from the edge
        newVertex->newPosition = e->newPosition;
      }
      e = nextEdge;
    }
    
    // 4. Flip any new edge that connects an old and new vertex.
    for (EdgeIter activeEdge = mesh.edgesBegin(); activeEdge != mesh.edgesEnd(); activeEdge++) {
      VertexIter v0 = activeEdge->halfedge()->vertex();
      VertexIter v1 = activeEdge->halfedge()->twin()->vertex();

      if (((v0->isNew && !v1->isNew) || (!v0->isNew && v1->isNew)) && activeEdge->isNew) {       
        mesh.flipEdge(activeEdge);
      }
    }

    // 5. Copy the new vertex positions into final Vertex::position.
    for (VertexIter activeVertex = mesh.verticesBegin(); activeVertex != mesh.verticesEnd(); activeVertex++) {
      activeVertex->position = activeVertex->newPosition;
    }
  }
}
