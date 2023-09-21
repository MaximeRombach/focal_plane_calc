
using SolidWorks.Interop.sldworks;

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

    }
}
