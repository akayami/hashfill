package hashfill

import (
	"github.com/mmcloughlin/geohash"
	"github.com/paulmach/orb"
	"github.com/paulmach/orb/planar"
	geom "github.com/twpayne/go-geom"
)

// Container tests if a hash is contained.
type Container interface {
	Contains(*geom.Polygon, string) (bool, error)
}

// Intersector tests if a hash intersects.
type Intersector interface {
	Intersects(*geom.Polygon, string) (bool, error)
}

type predicate func(geofence *geom.Polygon, hash string) (bool, error)

type containsFunc predicate

func (f containsFunc) Contains(p *geom.Polygon, hash string) (bool, error) {
	return f(p, hash)
}

type intersectsFunc predicate

func (f intersectsFunc) Intersects(p *geom.Polygon, hash string) (bool, error) {
	return f(p, hash)
}

// Intersects tests if the geofence contains the hash by doing a geos intersection.
var Intersects = intersectsFunc(func(geofence *geom.Polygon, hash string) (bool, error) {
	hashGeo := hashToGeometry(hash)
	fence := ConvertGeomPolygonToOrbPolygon(geofence)
	return polygonIntersectsPolygon(fence, hashGeo), nil
})

// Contains tests if the geofence contains the hash by doing a geos contains.
var Contains = containsFunc(func(geofence *geom.Polygon, hash string) (bool, error) {
	hashGeo := hashToGeometry(hash)
	fence := ConvertGeomPolygonToOrbPolygon(geofence)
	return polygonContainsPolygon(hashGeo, fence), nil
})

//func geomToGeosCoord(coord geom.Coord) geos.Coord {
//	return geos.Coord{
//		X: coord.X(),
//		Y: coord.Y(),
//	}
//}

//func geomToGeosCoords(coords []geom.Coord) []geos.Coord {
//	out := make([]geos.Coord, len(coords))
//	for i := 0; i < len(coords); i++ {
//		out[i] = geomToGeosCoord(coords[i])
//	}
//	return out
//}

// hashToGeometry converts a a geohash to a geos polygon by taking its bounding box.
func hashToGeometry(hash string) orb.Polygon {
	bounds := geohash.BoundingBox(hash)
	return orb.Polygon{
		{
			{bounds.MinLng, bounds.MinLat},
			{bounds.MinLng, bounds.MaxLat},
			{bounds.MaxLng, bounds.MaxLat},
			{bounds.MaxLng, bounds.MinLat},
			{bounds.MinLng, bounds.MinLat},
		},
	}
}

//func polygonToGeometry(geofence *geom.Polygon) *geos.Geometry {
//	// Convert the outer shell to geos format.
//	shell := geofence.LinearRing(0).Coords()
//	shellGeos := geomToGeosCoords(shell)
//
//	// Convert each hole to geos format.
//	numHoles := geofence.NumLinearRings() - 1
//	holes := make([][]geos.Coord, numHoles)
//	for i := 0; i < numHoles; i++ {
//		holes[i] = geomToGeosCoords(geofence.LinearRing(i).Coords())
//	}
//
//	return geos.Must(geos.NewPolygon(shellGeos, holes...))
//}

// ConvertGeomPolygonToOrbPolygon converts a geom.Polygon to an orb.Polygon
func ConvertGeomPolygonToOrbPolygon(gPolygon *geom.Polygon) orb.Polygon {
	numRings := gPolygon.NumLinearRings()
	oPolygon := make(orb.Polygon, numRings)

	for i := 0; i < numRings; i++ {
		ring := gPolygon.LinearRing(i)
		coords := ring.FlatCoords()
		numCoords := len(coords) / 2 // Divide by 2 as each coordinate is a pair of values
		oRing := make(orb.Ring, numCoords)

		for j := 0; j < numCoords; j++ {
			oRing[j] = orb.Point{coords[2*j], coords[2*j+1]}
		}

		oPolygon[i] = oRing
	}

	return oPolygon
}

func polygonContainsPolygon(outer, inner orb.Polygon) bool {
	for _, ring := range inner {
		for _, point := range ring {
			if !planar.PolygonContains(outer, point) {
				return false
			}
		}
	}
	return true
}

func polygonIntersectsPolygon(outer, inner orb.Polygon) bool {
	if !inner.Bound().Intersects(outer.Bound()) {
		return false
	}
	if !outer.Bound().Intersects(inner.Bound()) {
		return false
	}

	return polygonsIntersect(outer, inner)
}

func lineIntersects(a1, a2, b1, b2 orb.Point) bool {
	// Based on the "Line segment intersection" algorithm
	det := (a2[0]-a1[0])*(b2[1]-b1[1]) - (b2[0]-b1[0])*(a2[1]-a1[1])
	if det == 0 {
		return false // Collinear or parallel
	}
	lambda := ((b2[1]-b1[1])*(b2[0]-a1[0]) + (b1[0]-b2[0])*(b2[1]-a1[1])) / det
	gamma := ((a1[1]-a2[1])*(b2[0]-a1[0]) + (a2[0]-a1[0])*(b2[1]-a1[1])) / det
	return (0 < lambda && lambda < 1) && (0 < gamma && gamma < 1)
}

func polygonsIntersect(poly1, poly2 orb.Polygon) bool {
	for _, ring1 := range poly1 {
		for i := 0; i < len(ring1)-1; i++ {
			edge1Start, edge1End := ring1[i], ring1[i+1]
			for _, ring2 := range poly2 {
				for j := 0; j < len(ring2)-1; j++ {
					edge2Start, edge2End := ring2[j], ring2[j+1]
					if lineIntersects(edge1Start, edge1End, edge2Start, edge2End) {
						return true
					}
				}
			}
		}
	}
	return false
}
