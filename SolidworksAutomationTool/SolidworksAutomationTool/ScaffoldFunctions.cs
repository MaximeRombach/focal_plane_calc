
using SolidWorks.Interop.sldworks;
using SolidWorks.Interop.swconst;
using System.Diagnostics;

namespace SolidworksAutomationTool
{
    public class ScaffoldFunctions
    {
        /* Display a prompt to the console and wait for the user's input before continuing */
        public static void PromptAndWait(string prompt)
        {
            Console.WriteLine(prompt);
            Console.ReadKey();
        }
        public static double DegreeToRadian(double angleInDegrees)
        {
            return angleInDegrees / 180.0 * Math.PI;
        }

        // some scaffold functions. So far none of the functions clears selection. So the user should remember to manually clear selection
        /*
            Wrapper function to make selected segments coincident
            Note that this function doesn't check if the selected segments are valid to be make coincident with each others.
         */
        public static void MakeSelectedCoincide(ref ModelDoc2 partModelDoc)
            => partModelDoc.SketchAddConstraints("sgCOINCIDENT");

        /*
            Wrapper function to make a selected line vertical. 
            Note that this function doesn't check if the selected is a line or something else. 
            TODO: add check to make sure the selected item is a valid line
         */
        public static void MakeSelectedLineVertical(ref ModelDoc2 partModelDoc)
            => partModelDoc.SketchAddConstraints("sgVERTICAL2D");

        /*
            Wrapper function to make a selected line horizontal. 
            Note that this function doesn't check if the selected is a line or something else. 
            TODO: add check to make sure the selected item is a valid line
         */
        public static void MakeSelectedLineHorizontal(ref ModelDoc2 partModelDoc)
            => partModelDoc.SketchAddConstraints("sgHORIZONTAL2D");

        // Wrapper function to clear selection
        public static void ClearSelection(ref ModelDoc2 partModelDoc) 
            => partModelDoc.ClearSelection2(true);

        // Wrapper function to select the origin
        public static void SelectOrigin(ref ModelDoc2 partModelDoc)
            => partModelDoc.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, false, 0, null, 0);

        /* Wrapper function to zoom-to-fit the view */
        public static void ZoomToFit(ref ModelDoc2 partModelDoc)
            => partModelDoc.ViewZoomtofit2();

        /* A debug function that traverses the sketch segments in a polygon and prints the type of each sketch segment
         * There's no built-in funciton to get see the structure of the polygon. 
         */
        public static void PrintPolygonDataStructure(ref object[] polygon)
        {
            Debug.WriteLine($"Found {polygon.Length} sketch segments in the polygon");
            foreach (SketchSegment triangleSegment in polygon)
            {
                string triangleSegmentType = triangleSegment.GetType() switch
                {
                    (int)swSketchSegments_e.swSketchARC => triangleSegmentType = "Arc",
                    (int)swSketchSegments_e.swSketchLINE => triangleSegmentType = "Line",
                    (int)swSketchSegments_e.swSketchELLIPSE => triangleSegmentType = "Ellipse",
                    (int)swSketchSegments_e.swSketchPARABOLA => triangleSegmentType = "Parabola",
                    (int)swSketchSegments_e.swSketchSPLINE => triangleSegmentType = "Spline",
                    // Actually, why a sketch segment can be of type text??
                    (int)swSketchSegments_e.swSketchTEXT => triangleSegmentType = "Text",
                    _ => "Unknown",
                };
                Debug.WriteLine($"Segment is of type: {triangleSegmentType}");
            }
        }

        /* A function to get the center point of the inscribed construction circle inside the triangle polygon
         * Returns the center point as a sketch point if the polygon contains a Sketch Arc
         *          else, returns null.
         */
        public static SketchPoint? GetTriangleCenterPoint(ref object[] polygon)
        {
            foreach (SketchSegment triangleSegment in polygon.Reverse().Cast<SketchSegment>())
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchARC)
                {
                    return (SketchPoint)((SketchArc)triangleSegment).GetCenterPoint2(); ;
                }
            }
            return null;
        }

        /* A function to get one of the sides of a triangle polygon
         * Returns a side of the triangle if the polygon contains at least a Sketch line
         *  else Returns null
         */
        public static SketchLine? GetOneTriangleSide(ref object[] polygon)
        {
            foreach(SketchSegment triangleSegment in polygon.Cast<SketchSegment>())
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    return (SketchLine)triangleSegment;
                }
            }
            return null;
        }

        /* A wrapper function to reduce the boilerplate code for creating normal planes using the "point and normal" method
         *  This function first selects an extrusion axis, requiring it to be perpendicular to the ref plane;
         *  then selects a point which provides the "offset" of the plane, requiring it to be coincident with the ref plane
         *  NOTE: this method does NOT change the selection list. Users should manually clear the selections.
         */
        public static RefPlane CreateRefPlaneFromPointAndNormal(SketchPoint point, SketchSegment normal, SelectData swSelectData, FeatureManager featureManager)
        {
            swSelectData.Mark = 0;
            normal.Select4(true, swSelectData);
            swSelectData.Mark = 1;
            point.Select4(true, swSelectData);
            // although the InsertRefPlane function can take 3 constraints, we don't need to provide the third constraint,
            // as the 2 constraints are enough for the "point and normal" method
            RefPlane refPlane = (RefPlane)featureManager.InsertRefPlane((int)swRefPlaneReferenceConstraints_e.swRefPlaneReferenceConstraint_Perpendicular, 0,
                                                                                          (int)swRefPlaneReferenceConstraints_e.swRefPlaneReferenceConstraint_Coincident, 0, 0, 0);
            return refPlane;
        }

        /* Wrapper function to add dimension to the selected object. 
         * Returns the model dimension, NOT the DisplayDimension
         * NOTE: this function does NOT change the selection list. Users should manually clear the selections.
         * NOTE: this function can NOT be used to dimension chamfers.
         */
        public static Dimension AddDimensionToSelected(ref ModelDoc2 partModelDoc, double dimensionValue, double dimLocationX, double dimLocationY, double dimLocationZ)
        {
            // first create a DisplayDimension that 
            DisplayDimension displayDimension = (DisplayDimension)partModelDoc.AddDimension2(dimLocationX, dimLocationY, dimLocationZ);
            // The Index argument is valid for chamfer display dimensions only.
            // If the display dimension is not a chamfer display dimension, then Index is ignored.
            // To get both chamfer display dimensions, you must call this property twice; specify 0 for Index in the first call and 1 for Index in the second call.
            Dimension dimension = displayDimension.GetDimension2(0);
            // not sure if swSetValue_InThisConfiguration is the best parameter to use
            dimension.SetSystemValue3(dimensionValue, (int)swSetValueInConfiguration_e.swSetValue_InThisConfiguration, "");
            return dimension;
        }

        /* Overloaded AddDimensionToSelected function to simplify function calls
         * Often times we already have a sketch point indicating a position. We should be able to simply pass the point in as a position indicator, instead of having to specify x,y,z coordinates manually
         */
        public static Dimension AddDimensionToSelected(ref ModelDoc2 partModelDoc, double dimensionValue, SketchPoint displayDimensionLocation )
        {
            // first create a DisplayDimension that 
            DisplayDimension displayDimension = (DisplayDimension)partModelDoc.AddDimension2(displayDimensionLocation.X, displayDimensionLocation.Y, displayDimensionLocation.Z);
            // The Index argument is valid for chamfer display dimensions only.
            // If the display dimension is not a chamfer display dimension, then Index is ignored.
            // To get both chamfer display dimensions, you must call this property twice; specify 0 for Index in the first call and 1 for Index in the second call.
            Dimension dimension = displayDimension.GetDimension2(0);
            // not sure if swSetValue_InThisConfiguration is the best parameter to use
            dimension.SetSystemValue3(dimensionValue, (int)swSetValueInConfiguration_e.swSetValue_InThisConfiguration, "");
            return dimension;
        }

        /* A function to make a block of chamfered triangle from a triangle polygon
         */
        public static SketchBlockDefinition MakeChamferedTriangleBlockFromTrianglePolygon(object[] trianglePolygon, double chamferLength, ref ModelDoc2 partModelDoc, SelectData swSelectData)
        {
            // get the vertices of the triangle. The hashSet will pick only the unique vertices
            HashSet<SketchPoint> verticesInTriangleSet = new();
            foreach (SketchSegment triangleSegment in trianglePolygon.Cast<SketchSegment>())
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    verticesInTriangleSet.Add((SketchPoint)((SketchLine)triangleSegment).GetStartPoint2());
                    verticesInTriangleSet.Add((SketchPoint)((SketchLine)triangleSegment).GetEndPoint2());
                }
            }
            ClearSelection(ref partModelDoc);
            // make a chamfer at every vertex
            List<SketchSegment> chamferSegments = new(3);
            foreach (SketchPoint vertex in verticesInTriangleSet)
            {
                // make chamfers
                vertex.Select4(true, swSelectData);
                SketchSegment chamferSegment = partModelDoc.SketchManager.CreateChamfer((int)swSketchChamferType_e.swSketchChamfer_DistanceEqual, chamferLength, chamferLength);
                chamferSegments.Add(chamferSegment);
            }
            // select all chamfer segments
            chamferSegments.ForEach(chamferSegment => chamferSegment.Select4(true, swSelectData));

            // Select the triangle polygon
            foreach (SketchSegment triangleSegment in trianglePolygon.Cast<SketchSegment>())
            {
                triangleSegment.Select4(true, swSelectData);
            }

            // select the 3 vertices of the original triangle
            verticesInTriangleSet.ToList<SketchPoint>().ForEach(triangleVertex => triangleVertex.Select4(true, swSelectData));

            // Great! The block was sucessfully created. NOTE: for some reasons, MultiSelect2 Method (IModelDocExtension) does NOT select anything in the triangle polygon
            SketchBlockDefinition chamferedTriangleBlock = partModelDoc.SketchManager.MakeSketchBlockFromSelected(null);

            return chamferedTriangleBlock;
        }

        /* A function to make a block of triangle from a triangle polygon
         */ 
        public static SketchBlockDefinition MakeTriangleBlockFromTrianglePolygon(object[] trianglePolygon, ref ModelDoc2 partModelDoc, SelectData swSelectData)
        {
            // Select the triangle polygon
            foreach (SketchSegment triangleSegment in trianglePolygon.Cast<SketchSegment>())
            {
                triangleSegment.Select4(true, swSelectData);
            }
            SketchBlockDefinition triangleBlock = partModelDoc.SketchManager.MakeSketchBlockFromSelected(null);
            return triangleBlock;
        }
    }
}