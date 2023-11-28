
using SolidWorks.Interop.sldworks;
using SolidWorks.Interop.swconst;
using System.Diagnostics;

namespace SolidworksAutomationTool
{
    /* A struct storing the 3 basic reference planes */
    public struct BasicReferenceGeometry
    {
        public Feature frontPlane;
        public Feature topPlane;
        public Feature rightPlane;
        public SketchPoint origin;
    }

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

        /* Wrapper function to make selected two lines parallel - Untested
         * NOTE: this function does NOT check if the selected objects are just two lines or something else
         * NOTE: this function does not clear the selection. The user should manually clear selection
         */
        public static void MakeSelectedLinesParallel(ref ModelDoc2 partModelDoc)
            => partModelDoc.SketchAddConstraints("sgPARALLEL");

        // Wrapper function to clear selection
        public static void ClearSelection(ref ModelDoc2 partModelDoc)
            => partModelDoc.ClearSelection2(true);

        // Wrapper function to select the origin. ONLY WORKS IN THE ENGLISH VERSION OF SOLIDWORKS
        public static void SelectOrigin(ref ModelDoc2 partModelDoc)
            => partModelDoc.Extension.SelectByID2("Point1@Origin", "EXTSKETCHPOINT", 0, 0, 0, false, 0, null, 0);

        /* Wrapper function to select the sketch with the given name */
        public static void SelectSketch(ref ModelDoc2 partModelDoc, string sketchName, bool appendToSelection = false)
        {
            partModelDoc.Extension.SelectByID2(sketchName, "SKETCH", 0, 0, 0, appendToSelection, 0, null, 0);
        }

        /* Wrapper function to select all sketch segments in a given array.
         * NOTE: this function does not clear previous selections. The user should manually clear selections if needed
         */
        public static void SelectAllSketchSegments(ref object[] segmentArray, SelectData swSelectData)
        {
            foreach (SketchSegment segment in segmentArray.Cast<SketchSegment>())
            {
                segment.Select4(true, swSelectData);
            }
        }

        /*Wrapper function to rotate selected sketch segments by certain angle.
         * NOTE: this function preserves the sketch relations from the source sketch
         */
        public static void RotateSelected(ref ModelDoc2 partModelDoc, double rotationCenterX, double rotationCenterY, double angleInRad)
        {
            // preserve the sketch relations when copying a sketch
            partModelDoc.Extension.RotateOrCopy(false, 1, true, rotationCenterX, rotationCenterY, 0, 0, 0, 1, angleInRad);
        }

        /* Wrapper function to zoom-to-fit the view */
        public static void ZoomToFit(ref ModelDoc2 partModelDoc)
            => partModelDoc.ViewZoomtofit2();

        /* A debug function that traverses the sketch segments in a polygon and prints the type of each sketch segment
         * There's no built-in funciton to get see the structure of the polygon. 
         */
        public static void PrintPolygonDataStructure(ref object[] polygon)
        {
            Debug.WriteLine($"Found {polygon.Length} sketch segments in the polygon");
            foreach (SketchSegment triangleSegment in polygon.Cast<SketchSegment>())
            {
                string triangleSegmentType = triangleSegment.GetType() switch
                {
                    (int)swSketchSegments_e.swSketchARC => "Arc",
                    (int)swSketchSegments_e.swSketchLINE => "Line",
                    (int)swSketchSegments_e.swSketchELLIPSE => "Ellipse",
                    (int)swSketchSegments_e.swSketchPARABOLA => "Parabola",
                    (int)swSketchSegments_e.swSketchSPLINE => "Spline",
                    // Actually, why a sketch segment can be of type text??
                    (int)swSketchSegments_e.swSketchTEXT => "Text",
                    _ => "Unknown",
                };
                Debug.WriteLine($"Segment is of type: {triangleSegmentType}");
            }
        }

        /* Print a sketch point's coordinates to the Debug stream */
        public static void PrintSketchPoint(SketchPoint aPoint, string pointName)
        {
            Debug.WriteLine($"Point {pointName}: x: {aPoint.X}, y: {aPoint.Y}, z: {aPoint.Z}");
        }

        /* Print a math point's coordinates to the Debug stream */
        public static void PrintMathPoint(MathPoint aPoint, string pointName)
        {
            double[] pointDataArray = (double[])aPoint.ArrayData;
            Debug.WriteLine($"Point {pointName}: x: {pointDataArray[0]}, y: {pointDataArray[1]}, z: {pointDataArray[2]}");
        }

        /* Debug function to see the features inside the feature manager design tree
     * This function can be used to check if a feature is created as expected
     */
        public static void PrintFeaturesInFeatureManagerDesignTree(ref ModelDoc2 partModelDoc)
        {
            Debug.WriteLine("Printing features in this part:");
            int numberOfFeatures = partModelDoc.GetFeatureCount();
            Debug.WriteLine($"Found {numberOfFeatures} features in this design tree");
            // Solidworks doesn't seem to have an api to traverse the design tree from the beginning. So we do some reverse indexing here
            for (int featureIdx = numberOfFeatures - 1; featureIdx >= 0; featureIdx--)
            {
                Feature aFeature = (Feature)partModelDoc.FeatureByPositionReverse(featureIdx);
                // format the print a bit to make it easy to read
                Debug.WriteLine($"feature #{numberOfFeatures -1 - featureIdx, 5}    type: {aFeature.GetTypeName2(), 25},    name: {aFeature.Name, 25}");
            }
        }

        /* A wrapper function to get the name of the current active sketch*/
        public static string GetActiveSketchName(ref ModelDoc2 partModelDoc)
        {
            // TODO: add a check on the active status of the active sketch. There could be no active sketch
            return ((Feature)partModelDoc.SketchManager.ActiveSketch).Name;
        }

        /* Get the under lying basic reference geometry: front, top, right, and the origin.
         * This function tries to get the basic reference planes without using their names to avoid unexpected behavior on computers using other languages than English
         * For example, when selecting planes using their english names, french-based solidworks will fail to locate the planes
         */
        public static BasicReferenceGeometry GetBasicReferenceGeometry(ref ModelDoc2 partModelDoc)
        {
            BasicReferenceGeometry basicReferenceGeometry = new();
            // initalize the reference geometry counters to keep track of the number of found ref geometry. 
            uint refPlanesLocated = 0;
            uint originLocated = 0;
            // Solidworks doesn't seem to have an api to traverse the design tree by index from the beginning. So we do some reverse indexing here
            for (int featureIdx = partModelDoc.GetFeatureCount() - 1; featureIdx >= 0; featureIdx--)
            {
                // no need to continue looking into features if already found 3 basic ref planes and 1 origin
                if (refPlanesLocated >= 3 && originLocated >= 1)
                {
                    break;
                }
                // check if the current feature is of type RefPlane. If so, update the basicReferencePlanes struct
                Feature aFeature = (Feature)partModelDoc.FeatureByPositionReverse(featureIdx);
                if (aFeature.GetTypeName2() == "RefPlane")
                {
                    switch (refPlanesLocated)
                    {
                        case 0:
                            basicReferenceGeometry.frontPlane = aFeature;
                            refPlanesLocated += 1;
                            Debug.WriteLine($"default front plane name: {basicReferenceGeometry.frontPlane?.Name}");
                            break;
                        case 1:
                            basicReferenceGeometry.topPlane = aFeature;
                            refPlanesLocated += 1;
                            Debug.WriteLine($"default top plane name: {basicReferenceGeometry.topPlane?.Name}");
                            break;
                        case 2:
                            basicReferenceGeometry.rightPlane = aFeature;
                            refPlanesLocated += 1;
                            Debug.WriteLine($"default right plane name: {basicReferenceGeometry.rightPlane?.Name}");
                            break;
                    }
                }
                // origin has the type name "OriginProfileFeature", at least in solidworks 2023
                else if (aFeature.GetTypeName2() == "OriginProfileFeature")
                {
                    if (originLocated == 0)
                    {
                        Debug.WriteLine($"default origin name: {aFeature.Name}");
                        // the surprise here is that the origin is actually a sketch point, instead of a feature. So we need to get the sketch point underneath the "origin feature"
                        object[] listOfSketchPoints = (object[])((Sketch)aFeature.GetSpecificFeature2()).GetSketchPoints2();
                        basicReferenceGeometry.origin = (SketchPoint)listOfSketchPoints[0];
                        originLocated += 1;
                    }
                }
            }
            return basicReferenceGeometry;
        }

        /* Get the distance between two math points */
        public static double GetDistanceBetweenTwoMathPoints(MathPoint p1, MathPoint p2)
        {
            double[] p1DataArray = (double[])p1.ArrayData;
            double[] p2DataArray = (double[])p2.ArrayData;

            return Math.Sqrt(Math.Pow(p1DataArray[0] - p2DataArray[0], 2.0) +
                                Math.Pow(p1DataArray[1] - p2DataArray[1], 2.0) +
                                Math.Pow(p1DataArray[2] - p2DataArray[2], 2.0));
        }

        /* Get the distance between two sketch points */
        public static double GetDistanceBetweenTwoSketchPoints(SketchPoint p1, SketchPoint p2)
        {
            return Math.Sqrt(   Math.Pow(p1.X - p2.X, 2.0) +
                                Math.Pow(p1.Y - p2.Y, 2.0) +
                                Math.Pow(p1.Z - p2.Z, 2.0));
        }

        /* Get the number of features in this document */
        public static int GetFeatureCount(ref ModelDoc2 partModelDoc)
            => partModelDoc.FeatureManager.GetFeatureCount(false);

        /* Get the index of the closest sketch point to the origin. 
         * Params: sketchPoints: reference to a list of sketchPoints
         * Returns: the index of the closest sketch point
         * This function looks for the closest sketchpoint based on the shortest Euclidean distance in 3D (aka L2 norm of 3D vector) from the origin
         */
        public static int GetIndexSketchPointClosestToOrigin(ref List<SketchPoint> sketchPoints)
        {
            double shortestDistance = double.MaxValue;
            int closestPointIdx = -1;
            for (int idx = 0; idx < sketchPoints.Count; idx++)
            {
                SketchPoint point = sketchPoints[idx];
                // not taking the square root to save computation cycles
                double distance = (point.X * point.X + point.Y * point.Y + point.Z * point.Z);
                if (distance < shortestDistance)
                {
                    shortestDistance = distance;
                    closestPointIdx = idx;
                }
            }
            return closestPointIdx;
        }

        /* A function to get the center point of the inscribed construction circle inside the triangle polygon
         * Returns the center point as a sketch point if the polygon contains a Sketch Arc
         *          else, returns null.
         */
        public static SketchPoint? GetTriangleCenterPoint(ref object[] trianglePolygon)
        {

            // the trick is to check the ForConstruction property of the circles. Only the inscribed circle is marked for construction
            // the pin holes are all solid shapes
            foreach (SketchSegment triangleSegment in trianglePolygon.Cast<SketchSegment>())
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchARC && triangleSegment.ConstructionGeometry)
                {
                    return (SketchPoint)((SketchArc)triangleSegment).GetCenterPoint2(); ;
                }
            }
            return null;
        }

        /* Returns a list of the vertices in a full triangle */
        public static List<SketchPoint> GetVerticesInTriangle(ref object[] trianglePolygon)
        {
            // get the vertices of the triangle.
            // the same vertex will be shared by two lines. We use the hashSet to pick only the unique vertices
            HashSet<SketchPoint> verticesInTriangleSet = new();
            foreach (SketchSegment triangleSegment in trianglePolygon.Cast<SketchSegment>())
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    verticesInTriangleSet.Add((SketchPoint)((SketchLine)triangleSegment).GetStartPoint2());
                    verticesInTriangleSet.Add((SketchPoint)((SketchLine)triangleSegment).GetEndPoint2());
                }
            }
            return verticesInTriangleSet.ToList();
        }

        /* A function to get one of the sides of a triangle polygon
         * Returns a side of the triangle if the polygon contains at least a Sketch line
         *  else Returns null
         */
        public static SketchLine? GetOneTriangleSide(ref object[] polygon)
        {
            foreach (SketchSegment triangleSegment in polygon.Cast<SketchSegment>())
            {
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    return (SketchLine)triangleSegment;
                }
            }
            return null;
        }

        /* Get the most horizontal side of a triangle polygon
         * TODO: add more descriptions on the trick used
         * */
        public static SketchLine? GetMostHorizontalTriangleSide(ref object[] polygon)
        {
            SketchLine? mostHorizontalSide = null;
            double smallestSlopeMagnitude = double.MaxValue;
            foreach (SketchSegment triangleSegment in polygon.Cast<SketchSegment>())
            {
                // calculate the slope of the projection of the side on the xy plane.
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    SketchPoint startPoint = (SketchPoint)((SketchLine)triangleSegment).GetStartPoint2();
                    SketchPoint endPoint = (SketchPoint)((SketchLine)triangleSegment).GetEndPoint2();
                    double sideSlopeMagnitude = Math.Abs((endPoint.Y - startPoint.Y) / (endPoint.X - startPoint.X));
                    if (sideSlopeMagnitude < smallestSlopeMagnitude)
                    {
                        smallestSlopeMagnitude = sideSlopeMagnitude;
                        mostHorizontalSide = (SketchLine)triangleSegment;
                    }
                }
            }
            return mostHorizontalSide;
        }

        /* Get the longest horizontal side of a chamfered triangle polygon
         */
        public static SketchLine? GetLongestMostHorizontalTriangleSide(ref object[] chamferedTrianglePolygon)
        {
            SketchLine? mostHorizontalLongSide = null;
            double smallestSlopeMagnitude = double.MaxValue;
            double longestSideLength = -1;
            foreach (SketchSegment side in chamferedTrianglePolygon.Cast<SketchSegment>())
            {
                // calculate the slope of the projection of the side on the xy plane.
                if (side.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    SketchPoint startPoint = (SketchPoint)((SketchLine)side).GetStartPoint2();
                    SketchPoint endPoint = (SketchPoint)((SketchLine)side).GetEndPoint2();
                    double sideSlopeMagnitude = Math.Abs((endPoint.Y - startPoint.Y) / (endPoint.X - startPoint.X));

                    double currentSideLength = side.GetLength();
                    // the slope of the current side is not small enough, continue to check other sides
                    if (sideSlopeMagnitude > smallestSlopeMagnitude)
                    {
                        continue;
                    }
                    // found a flatter side than all previous sides. Update the flattest side, longeset side length and slope magnitude trackers
                    if (sideSlopeMagnitude < smallestSlopeMagnitude)
                    {
                        smallestSlopeMagnitude = sideSlopeMagnitude;
                        mostHorizontalLongSide = (SketchLine)side;
                    }
                    // found a side as flat as the flattest side found previously, check if it is longer than the previously found side
                    else if (sideSlopeMagnitude == smallestSlopeMagnitude && currentSideLength > longestSideLength)
                    {
                        mostHorizontalLongSide = (SketchLine)side;
                    }
                }
            }
            return mostHorizontalLongSide;
        }

        /* Get the side with the most positive slope (only considering the projection on the local skecth XY plane) in a chamfered triangle
         * This function is often used for finding the side used to setting parallel constraints to the first triangle in a pizza slice
         * This function is designed to be called when the chamfered triangle is relatively upright
         * Since the chamfered triangle has two side in parallel within itself, we will find two sides having similar slopes.
         * We choose the side with the most positive slope
         */
        public static SketchSegment? GetMostPositiveSlopedSideInChamferedTriangle(ref object[] polygon)
        {
            SketchSegment? mostPositiveSlopedSide = null;
            double mostPositiveSlope = -999;
            foreach (SketchSegment triangleSegment in polygon.Cast<SketchSegment>())
            {
                // calculate the slop of the projection of the side on the xy plane.
                if (triangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    SketchPoint startPoint = (SketchPoint)((SketchLine)triangleSegment).GetStartPoint2();
                    SketchPoint endPoint = (SketchPoint)((SketchLine)triangleSegment).GetEndPoint2();
                    double sideSlope = (endPoint.Y - startPoint.Y) / (endPoint.X - startPoint.X);
                    // assuming there will be a side with positive slope
                    if (sideSlope < 0 || sideSlope < mostPositiveSlope)
                    {
                        continue;
                    }
                    mostPositiveSlope = sideSlope;
                    mostPositiveSlopedSide = triangleSegment;
                }
            }
            return mostPositiveSlopedSide;
        }

        /* A function to get one of the longer sides of a chamfered triangle.
         * Returns a longer side of the chamfered triangle if the polygon contains at least a sketch line
         *  else Returns null
         */
        public static SketchLine? GetOneChamferedTriangleLongSide(ref object[] polygon)
        {
            SketchSegment? longSide = null;
            double longestSideLength = 0;
            foreach (SketchSegment chamferedTriangleSegment in polygon.Cast<SketchSegment>())
            {
                if (chamferedTriangleSegment.GetType() == (int)swSketchSegments_e.swSketchLINE)
                {
                    if (longestSideLength < chamferedTriangleSegment.GetLength())
                    {
                        longestSideLength = chamferedTriangleSegment.GetLength();
                        longSide = chamferedTriangleSegment;
                    }
                }
            }
            return (SketchLine?)longSide;
        }

        /* 
         * Create global variables in the equation manager
         * Params:  partModelDoc: reference to the ModelDoc2 object
         *          expression: global variable assignment in string
         * Returns Index of the new equation if successfully added, -1 if error occured
         * If adding a global variable assignment that already exists, this method returns an error.
         * */
        public static int CreateGlobalVariableInAllConfigs(ref ModelDoc2 partModelDoc, string variableName, double variableValue)
        {
            EquationMgr equationManager = partModelDoc.GetEquationMgr();
            // syntax is referenced from https://help.solidworks.com/2022/English/api/sldworksapi/Add_Equations_Example_CSharp.htm?verRedirect=1
            // by default adding the global variable to the end of the equation list
            //int variableIndex = equationManager.Add3(-1, $"\"{variableName}\" = {variableValue}m", true, (int)swInConfigurationOpts_e.swAllConfiguration, null);
            string variableAssignmentExpression = $"\"{variableName}\" = {variableValue,0:F6}m";
            int variableIndex = equationManager.Add2(-1, variableAssignmentExpression, true);
            // TODO: check if the feature manager update is necessary
            partModelDoc.FeatureManager.UpdateFeatureTree();
            if (variableIndex == -1)
            {
                Console.WriteLine($"ERROR adding global variable : {variableName} with value {variableValue}");
            }
            return variableIndex;
        }

        /* A wrapper function to reduce the boilerplate code for creating normal planes using the "point and normal" method
         *  This function first selects an extrusion axis, requiring it to be perpendicular to the ref plane;
         *  then selects a point which provides the "offset" of the plane, requiring it to be coincident with the ref plane
         *  NOTE: this method does NOT change the selection list. Users should manually clear the selections.
         */
        public static RefPlane CreateRefPlaneFromPointAndNormal(SketchPoint point, SketchSegment normal, string? planeName, SelectData swSelectData, FeatureManager featureManager)
        {
            swSelectData.Mark = 0;
            normal.Select4(true, swSelectData);
            swSelectData.Mark = 1;
            point.Select4(true, swSelectData);
            // although the InsertRefPlane function can take 3 constraints, we don't need to provide the third constraint,
            // as the 2 constraints are enough for the "point and normal" method
            RefPlane refPlane = (RefPlane)featureManager.InsertRefPlane((int)swRefPlaneReferenceConstraints_e.swRefPlaneReferenceConstraint_Perpendicular, 0,
                                                                                          (int)swRefPlaneReferenceConstraints_e.swRefPlaneReferenceConstraint_Coincident, 0, 0, 0);

            // The user can decide the name of the plane by passing a string. If null is passed in, the plane's name is left to solidworks to decide
            if (planeName != null)
            {
                ((Feature)refPlane).Name = planeName;
            }
            return refPlane;
        }

        /* A wrapper function to make two-way extrusion on a already selected sketch 
         * Reuturns the extrusion feature
         * NOTE: this function does not check if the selected is a valid sketch or something else.
         */
        public static Feature CreateTwoWayExtrusion(ref ModelDoc2 partModelDoc)
        {
            // the FeatureCut4 api takes a ton of arguments. This wrapper function is to simplify the calling process.
            // https://help.solidworks.com/2023/english/api/sldworksapi/solidworks.interop.sldworks~solidworks.interop.sldworks.ifeaturemanager~featurecut4.html?verRedirect=1
            Feature extrusionFeature = partModelDoc.FeatureManager.FeatureCut4(
                                    false,  // true for single ended cut, false for double-ended cut
                                    false,  // True to remove material outside of the profile of the flip side to cut, false to not
                                    false,  // True for Direction 1 to be opposite of the default direction
                                    (int)swEndConditions_e.swEndCondThroughAllBoth, // Termination type for the first end
                                    (int)swEndConditions_e.swEndCondThroughAllBoth, // Termination type for the second end 
                                    1, 1,   // depth of extrusion for 1st and 2nd end in meters
                                    false, false, // True allows a draft angle in the first/second direction, false does not allow drafting in the first/second direction
                                    false, false, // True for the first/second draft angle to be inward, false to be outward; only valid when Dchk1/Dchk2 is true
                                    1, 1,   // Draft angle for the first end; only valid when Dchk1 is true
                                    false, false, // If you chose to offset the first/second end condition from another face or plane, then true specifies offset in direction away from the sketch, false specifies offset from the face or plane in a direction toward the sketch
                                    false, false,
                                    false, true, true, true, true, false,
                                    (int)swStartConditions_e.swStartSketchPlane,  // Start conditions as defined in swStartConditions_e
                                    0,      // If T0 is swStartConditions_e.swStartOffset, then specify an offset value
                                    false,  // If T0 is swStartConditions_e.swStartOffset, then true to flip the direction of cut, false to not
                                    false);
            return extrusionFeature;
        }

        /* A wrapper function to make two-way extrusion on a already selected sketch. 
         * Direction 1 is extruded till the given point,
         * Direction 2 is extruded all through.
         */
        public static Feature CreateTwoWayExtrusionD1ToPointD2ThroughAll(ref ModelDoc2 partModelDoc, SketchPoint vertexD1ExtrudeTo, SelectData swSelectData)
        {
            // Magic number 32 is used when selecting up-to surface, up-to vertex, or offset-from surface.
            // Source: https://help.solidworks.com/2023/english/api/sldworksapi/SOLIDWORKS.Interop.sldworks~SOLIDWORKS.Interop.sldworks.IFeatureManager~FeatureRevolve2.html
            swSelectData.Mark = 32;
            vertexD1ExtrudeTo.Select4(true, swSelectData);
            // the FeatureCut4 api takes a ton of arguments. This wrapper function is to simplify the calling process.
            // https://help.solidworks.com/2023/english/api/sldworksapi/solidworks.interop.sldworks~solidworks.interop.sldworks.ifeaturemanager~featurecut4.html?verRedirect=1
            Feature extrusionFeature = partModelDoc.FeatureManager.FeatureCut4(
                                    false,  // true for single ended cut, false for double-ended cut
                                    false,  // True to remove material outside of the profile of the flip side to cut, false to not
                                    false,  // True for Direction 1 to be opposite of the default direction
                                    (int)swEndConditions_e.swEndCondUpToVertex,  // Termination type for the first end
                                    (int)swEndConditions_e.swEndCondThroughAll,     // Termination type for the second end 
                                    0.01, 1,   // depth of extrusion for 1st and 2nd end in meters
                                    false, false, // True allows a draft angle in the first/second direction, false does not allow drafting in the first/second direction
                                    false, false, // True for the first/second draft angle to be inward, false to be outward; only valid when Dchk1/Dchk2 is true
                                    1, 1,   // Draft angle for the first end; only valid when Dchk1 is true
                                    false, false, // If you chose to offset the first/second end condition from another face or plane, then true specifies offset in direction away from the sketch, false specifies offset from the face or plane in a direction toward the sketch
                                    false, false,
                                    false, true, true, true, true, false,
                                    (int)swStartConditions_e.swStartSketchPlane,  // Start conditions as defined in swStartConditions_e
                                    0,      // If T0 is swStartConditions_e.swStartOffset, then specify an offset value
                                    false,  // If T0 is swStartConditions_e.swStartOffset, then true to flip the direction of cut, false to not
                                    false);
            // should we reset the select mark?
            //swSelectData.Mark = -1;
            return extrusionFeature;
        }

        /* A wrapper function to make two-way extrusion on a already selected sketch. 
         * Direction 1 is extruded till the given distance,
         * Direction 2 is extruded all through.
         */
        public static Feature CreateTwoWayExtrusionD1ToDistanceD2ThroughAll(ref ModelDoc2 partModelDoc, double distanceD1ExtrudeTo)
        { // the FeatureCut4 api takes a ton of arguments. This wrapper function is to simplify the calling process.
            // https://help.solidworks.com/2023/english/api/sldworksapi/solidworks.interop.sldworks~solidworks.interop.sldworks.ifeaturemanager~featurecut4.html?verRedirect=1
            Feature extrusionFeature = partModelDoc.FeatureManager.FeatureCut4(
                                    false,  // true for single ended cut, false for double-ended cut
                                    false,  // True to remove material outside of the profile of the flip side to cut, false to not
                                    false,  // True for Direction 1 to be opposite of the default direction
                                    (int)swEndConditions_e.swEndCondBlind,  // Termination type for the first end
                                    (int)swEndConditions_e.swEndCondThroughAll,     // Termination type for the second end 
                                    distanceD1ExtrudeTo, distanceD1ExtrudeTo,   // depth of extrusion for 1st and 2nd end in meters
                                    false, false, // True allows a draft angle in the first/second direction, false does not allow drafting in the first/second direction
                                    false, false, // True for the first/second draft angle to be inward, false to be outward; only valid when Dchk1/Dchk2 is true
                                    1, 1,   // Draft angle for the first end; only valid when Dchk1 is true
                                    false, false, // If you chose to offset the first/second end condition from another face or plane, then true specifies offset in direction away from the sketch, false specifies offset from the face or plane in a direction toward the sketch
                                    false, false,
                                    false, true, true, true, true, false,
                                    (int)swStartConditions_e.swStartSketchPlane,  // Start conditions as defined in swStartConditions_e
                                    0,      // If T0 is swStartConditions_e.swStartOffset, then specify an offset value
                                    false,  // If T0 is swStartConditions_e.swStartOffset, then true to flip the direction of cut, false to not
                                    false);
            return extrusionFeature;
        }

        /* A function to copy and paste a reference sketch on the given plane and rotate if necessary
         * This function is just to save time and repetitive code
         */
        public static (Sketch, SketchPoint, object[]) CreateACopyAndRotateReferenceSketch(ref ModelDoc2 partModelDoc, RefPlane plane, SelectData swSelectData, string referenceSketchName, bool isUprightTriangle)
        {
            // copy and paste the reference sketch
            SelectSketch(ref partModelDoc, referenceSketchName);
            partModelDoc.EditCopy();
            ((Feature)plane).Select2(true, -1);
            partModelDoc.Paste();

            // manually update feature tree
            partModelDoc.FeatureManager.UpdateFeatureTree();
            // select the last pasted pin hole triangle sketch
            Feature pastedSketchFeature = (Feature)partModelDoc.FeatureByPositionReverse(0);
            pastedSketchFeature.Select2(false, -1);
            Sketch pastedSketch = (Sketch)pastedSketchFeature.GetSpecificFeature2();
            ClearSelection(ref partModelDoc);

            // edit the pasted pin hole sketch
            ((Feature)pastedSketch).Select2(false, -1);
            partModelDoc.EditSketch();

            // it's a bit weird to cast the sketch back and forth, but only the Feature object has the sketch's name
            //string pastedPinHoleTriangleSketchName = ((Feature)pastedSketch).Name;
            object[] sketchSegments = (object[])pastedSketch.GetSketchSegments();

            SketchPoint? triangleCenter = GetTriangleCenterPoint(ref sketchSegments);
            ClearSelection(ref partModelDoc);
            if (!isUprightTriangle)
            {
                // orientation flag is false, meaning the module should be upside-down
                SelectAllSketchSegments(ref sketchSegments, swSelectData);
                RotateSelected(ref partModelDoc, triangleCenter.X, triangleCenter.Y, Math.PI);
                ClearSelection(ref partModelDoc);
            }
            return (pastedSketch, triangleCenter, sketchSegments);
        }

        /* Create a copy of the pin hole sketch on a given plane using a reference pin hole sketch*/
        public static Sketch CreatePinHoleSketchFromReferenceSketch(ref ModelDoc2 partModelDoc, RefPlane plane, SelectData swSelectData, string pinHoleSketchName, bool isUprightTriangle, SketchPoint pointToCoincide, SketchSegment sideToParallel)
        {
            // copy and paste the pin hole triangle sketch
            // copy the full triangle sketch //
            (Sketch pastedPinHoleSketch, SketchPoint pinHoleTriangleCenter, object[] pinHoleTriangleSegments) = CreateACopyAndRotateReferenceSketch(
                                                                                                                        ref partModelDoc,
                                                                                                                        plane,
                                                                                                                        swSelectData,
                                                                                                                        pinHoleSketchName,
                                                                                                                        isUprightTriangle);
            // make the center of the pin hole triangle conincident with the extrusion axis
            pinHoleTriangleCenter.Select4(true, swSelectData);
            pointToCoincide.Select4(true, swSelectData);
            MakeSelectedCoincide(ref partModelDoc);
            ClearSelection(ref partModelDoc);

            // Make the flattest side parallel to the given reference side(should be the flattest side of a full triangle)
            SketchLine? aLongSidePinHoleTriangle = GetMostHorizontalTriangleSide(ref pinHoleTriangleSegments);
            ((SketchSegment)aLongSidePinHoleTriangle).Select4(true, swSelectData);
            sideToParallel.Select4(true, swSelectData);
            MakeSelectedLinesParallel(ref partModelDoc);
            ClearSelection(ref partModelDoc);

            // quit editing the pin hole triangle sketch
            // TODO: try not rebuilding now
            partModelDoc.InsertSketch2(false);
            ClearSelection(ref partModelDoc);
            return pastedPinHoleSketch;
        }

        /* Create a copy of the chamfered triangle sketch on a given plane using a reference chamfered triangle sketch */
        public static (Sketch, SketchPoint, SketchLine) CreateChamferedTriangleSketchFromReferenceSketch(ref ModelDoc2 partModelDoc, RefPlane plane, SelectData swSelectData, string chamferedTriangleSketchName, bool isUprightTriangle, SketchPoint pointToCoincide)
        {
            // copy the chamfered triangle sketch //
            (Sketch pastedChamferedSketch, SketchPoint chamferedTriangleCenterPoint, object[] chamferedTriangleSegments) = CreateACopyAndRotateReferenceSketch(
                                                                                                                        ref partModelDoc,
                                                                                                                        plane,
                                                                                                                        swSelectData,
                                                                                                                        chamferedTriangleSketchName,
                                                                                                                        isUprightTriangle);
            chamferedTriangleCenterPoint?.Select4(true, swSelectData);
            pointToCoincide.Select4(true, swSelectData);
            MakeSelectedCoincide(ref partModelDoc);
            ClearSelection(ref partModelDoc);

            // TESTING: make one of the sides to be in parallel to the most positively sloped side of the first reference chamfered triangle
            // TODO: check if can set one side of the chamfered triangle to be in parallel with the first chamfered triangle
            //SketchSegment? unchamferedTriangleMostPositiveSlopedSide = GetMostPositiveSlopedSideInChamferedTriangle(ref segments);
            // TODO: temporary workaround, simply set one of the sides to be in horizontal to fully define each chamfered triangle
            SketchLine aLongSideChamferedTriangle = GetLongestMostHorizontalTriangleSide(ref chamferedTriangleSegments);
            ((SketchSegment)aLongSideChamferedTriangle).Select4(true, swSelectData);
            MakeSelectedLineHorizontal(ref partModelDoc);
            ClearSelection(ref partModelDoc);

            // get a handle on the longest most horizontal side of a chamfered triangle.
            // This side will be used as a reference to set parallel constraint to a side of a full triangle
            //SketchLine? aLongSideChamferedTriangle = GetLongestMostHorizontalTriangleSide(ref chamferedTriangleSegments);

            // quit editing sketch
            // TODO: try not rebuilding now
            partModelDoc.InsertSketch2(false);

            return (pastedChamferedSketch, chamferedTriangleCenterPoint, aLongSideChamferedTriangle);
        }

        /* Create a copy of the full triangle sketch on a given plane using a reference full triangle sketch */
        public static (Sketch, SketchPoint, SketchSegment) CreateFullTriangleSketchFromReferenceSketch( ref ModelDoc2 partModelDoc, 
                                                                                                        RefPlane plane, 
                                                                                                        SelectData swSelectData, 
                                                                                                        string fullTriangleSketchName, 
                                                                                                        bool isUprightTriangle, 
                                                                                                        SketchPoint pointToCoincide, 
                                                                                                        SketchSegment sideToParallel)
        {
            // copy the full triangle sketch //
            (Sketch lastFullTriangleSketch, SketchPoint fullTriangleCenterPoint, object[] fullTriangleSegments) = CreateACopyAndRotateReferenceSketch(
                                                                                                                        ref partModelDoc, 
                                                                                                                        plane, 
                                                                                                                        swSelectData, 
                                                                                                                        fullTriangleSketchName, 
                                                                                                                        isUprightTriangle);

            fullTriangleCenterPoint.Select4(true, swSelectData);
            pointToCoincide.Select4(true, swSelectData);
            MakeSelectedCoincide(ref partModelDoc);
            ClearSelection(ref partModelDoc);

            // make the flattest side parallel to a flattest long side of the chamfered triangle
            SketchLine? aLongSideFullTriangle = GetMostHorizontalTriangleSide(ref fullTriangleSegments);
            sideToParallel.Select4(true, swSelectData);
            ((SketchSegment)aLongSideFullTriangle).Select4(true, swSelectData);
            MakeSelectedLinesParallel(ref partModelDoc);

            // quit editing the fully triangle sketch
            // TODO: try not rebuilding now
            partModelDoc.InsertSketch2(false);

            return (lastFullTriangleSketch, fullTriangleCenterPoint, (SketchSegment)aLongSideFullTriangle);
        }

        /* A wrapper function to enable the dimension dialog when an operation requires some input */
        public static void EnableInputDimensionByUser(ref SldWorks solidworks)
        {
            solidworks.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swInputDimValOnCreate, true);
        }

        /* A wrapper function to disable the dimension dialog when an operation requires user input */
        public static void DisableInputDimensionByUser(ref SldWorks solidworks)
        {
            solidworks.SetUserPreferenceToggle((int)swUserPreferenceToggle_e.swInputDimValOnCreate, false);
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

        /* A function to add chamfer to all vertices on a equilateral triangle.
         * This function add chamfers in the same way as the function MakeChamferedTriangleBlockFromTrianglePolygon but does not create a block
         * NOTE: this function WILL clear your previous selections
         * Return: void
         */
        public static void MakeChamferedTriangleFromTrianglePolygon(object[] trianglePolygon, double chamferLength, ref ModelDoc2 partModelDoc, SelectData swSelectData)
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
            foreach (SketchPoint vertex in verticesInTriangleSet)
            {
                // make chamfers
                vertex.Select4(true, swSelectData);
                SketchSegment chamferSegment = partModelDoc.SketchManager.CreateChamfer((int)swSketchChamferType_e.swSketchChamfer_DistanceEqual, chamferLength, chamferLength);
                ClearSelection(ref partModelDoc);
            }
        }

        /* A function to make a block of chamfered triangle from a triangle polygon. Currently not used
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
            verticesInTriangleSet.ToList().ForEach(triangleVertex => triangleVertex.Select4(true, swSelectData));

            // Great! The block was sucessfully created. NOTE: for some reasons, MultiSelect2 Method (IModelDocExtension) does NOT select anything in the triangle polygon
            SketchBlockDefinition chamferedTriangleBlock = partModelDoc.SketchManager.MakeSketchBlockFromSelected(null);

            return chamferedTriangleBlock;
        }

        /* A function to make a block of triangle from a triangle polygon. Currently not used
         */ 
        public static SketchBlockDefinition MakeTriangleBlockFromTrianglePolygon(object[] trianglePolygon, ref ModelDoc2 partModelDoc, SelectData swSelectData)
        {
            // Select the triangle polygon
            SelectAllSketchSegments(ref trianglePolygon, swSelectData);
            SketchBlockDefinition triangleBlock = partModelDoc.SketchManager.MakeSketchBlockFromSelected(null);
            return triangleBlock;
        }
    }
}