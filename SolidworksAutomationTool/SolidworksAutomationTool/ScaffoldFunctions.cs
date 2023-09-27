
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
        {
            partModelDoc.SketchAddConstraints("sgCOINCIDENT");
        }

        /*
            Wrapper function to make a selected line vertical. 
            Note that this function doesn't check if the selected is a line or something else. 
            TODO: add check to make sure the selected item is a valid line
         */
        public static void MakeSelectedLineVertical(ref ModelDoc2 partModelDoc)
        {
            partModelDoc.SketchAddConstraints("sgVERTICAL2D");
        }

        /*
            Wrapper function to make a selected line horizontal. 
            Note that this function doesn't check if the selected is a line or something else. 
            TODO: add check to make sure the selected item is a valid line
         */
        public static void MakeSelectedLineHorizontal(ref ModelDoc2 partModelDoc)
        {
            partModelDoc.SketchAddConstraints("sgHORIZONTAL2D");
        }

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
            foreach (SketchSegment triangleSegment in polygon.Reverse())
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
         */
        public static SketchLine? GetOneTriangleSide(ref object[] polygon)
        {
            foreach(SketchSegment triangleSegment in polygon)
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    return (SketchLine)triangleSegment;
                }
            }
            return null;
        }
    }
}
