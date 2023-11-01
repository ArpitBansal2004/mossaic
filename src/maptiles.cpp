/**
 * @file maptiles.cpp
 * Code for the maptiles function.
 */

#include <iostream>
#include <map>

#include "maptiles.h"

using namespace std;


Point<3> convertToXYZ(LUVAPixel pixel) {
    return Point<3>( pixel.l, pixel.u, pixel.v );
}

MosaicCanvas* mapTiles(SourceImage const& theSource,
                       vector<TileImage>& theTiles)
{
    map<Point<3>, TileImage*> conversion_map;
    vector<Point<3>> nodes;
    MosaicCanvas * canvas = new MosaicCanvas(theSource.getRows(), theSource.getColumns());
    
    for (auto it = theTiles.begin(); it != theTiles.end(); ++it) {
        Point<3> avg_color = convertToXYZ(it -> getAverageColor());
        conversion_map[avg_color] = &(*it);
        nodes.push_back(avg_color);
    }
    KDTree<3> kdtree(nodes);
    for (int row_tr = 0; row_tr < canvas->getRows(); ++row_tr) {
        for (int col_tr = 0; col_tr < canvas->getColumns(); ++col_tr) {
            Point<3> region_color = convertToXYZ(theSource.getRegionColor(row_tr, col_tr));
            Point<3> closest_color = kdtree.findNearestNeighbor(region_color);
            canvas->setTile(row_tr, col_tr, conversion_map[closest_color]);
        }
    }
    return canvas;
}

