
namespace SolidworksAutomationTool
{
    // this enum defines the allowed units in the txt file to be parsed in
    public enum Units: ushort
    {
        Millimeter = 0,
        Meter = 1,
    }
    /* Point Cloud class that supports reading generated a point cloud from a txt file */
    public class PointCloud
    {
        // All 3D points are stored in an array
        public List<Point3D> point3Ds = new List<Point3D>();

        /* The backbone of reading a point cloud from a txt file. 
         * Return:  false if the provided txt file does NOT exist / error occurred while reading the txt file.
         *          true if the txt file is read successfully. The resulting list of 3D points is stored in the class attribute point3Ds
         */

        // Simple points-printing function. Could be used for debugging
        public void PrintPoint3Ds()
        {
            // this number controls the field width (the fixed space taken in the console) for each number.
            const int FieldWidthRightAligned = 15;
            int pointIndex = -1;
            foreach (Point3D point in point3Ds)
            {
                pointIndex += 1;
                Console.WriteLine($"Point {pointIndex, 5}:      x: {point.x,FieldWidthRightAligned:F3}, y: {point.y,FieldWidthRightAligned:F3}, z: {point.z,FieldWidthRightAligned:F3}");
            }
        }

        /* 
         * Reads a txt file and parse all points into a list of 3D points
         * Returns: true if the parse was sucessful
         *          false otherwise
         */
        public bool ReadPointCloudFromTxt(string fileName, Units dataUnit)
        {
            // security check
            if ( !File.Exists(fileName) )
            {
                Console.WriteLine("ERROR: Cannot read point cloud from txt. Txt file doesn't exist");
                return false;
            }

            if ( !Enum.IsDefined( typeof( Units ), dataUnit ) )
            {
                Console.WriteLine("ERROR: the unit entered in <function>ReadPointCloudFromTxt is not valid");
                return false;
            }

            using (StreamReader streamReader = new StreamReader(fileName))
            {
                // line number starts from 0. 
                uint currentLine = 0;
                // define delimiters
                char[] splitOptions = { ' ' };

                // do a first read to go over the header. 
                string lineRead = streamReader.ReadLine();

                // start actually reading the first line of data
                lineRead = streamReader.ReadLine() ;
                currentLine += 1;

                while (lineRead != null)
                {
                    // split the line read by space
                    string[] splittedLine = lineRead.Split( splitOptions, StringSplitOptions.RemoveEmptyEntries );

                    // if a line contains anything other than 3 numbers, the data format is wrong. We need 3 points to define a point in 3D
                    if (splittedLine.Length != 3)
                    {
                        Console.WriteLine($"ERROR: Wrong data format in line {currentLine}. This line has {splittedLine.Length} numbers instead of 3");
                        return false;
                    }

                    // Try to convert the splitted line into 3 numbers and create a Point3D instance
                    float[] convertedNumbers = new float[3];
                    for (uint axis = 0; axis < splittedLine.Length; axis ++)
                    {
                        if (!Single.TryParse(splittedLine[axis], out convertedNumbers[axis]))
                        {
                            Console.WriteLine($"ERROR: Wrong data format in line {currentLine}. Please check for typos");
                            return false;
                        }

                        // Since solidworks's api seems to only take data in meters. We need to convert everything to meters, when the point cloud txt file contains data in mm
                        convertedNumbers[axis] = dataUnit switch
                        {
                            Units.Millimeter => convertedNumbers[axis] / 1000.0f,
                            Units.Meter => convertedNumbers[axis],
                            // Currently only support meters and milimeters. Other units might be added later
                            _ => -1,
                        };
                    }
                    // add the create point to the Point3D list
                    point3Ds.Add( new Point3D(convertedNumbers) );

                    // read the next line
                    lineRead = streamReader.ReadLine();
                    currentLine += 1;
                }
            }
            // in the case where the txt file is empty, let the user know.
            if (point3Ds.Count == 0)
            {
                Console.WriteLine("WARNING: no point cloud found in this txt file");
                return false;
            }
            Console.WriteLine($"Successfully read {point3Ds.Count} points from txt");
            return true;
        }
    }

    /* A simple class to hold a 3D point */
    public class Point3D : IEquatable<Point3D>
    {
        public float x = 0.0f;
        public float y = 0.0f;
        public float z = 0.0f;

        // Constructors of Point3D
        public Point3D() { }
        public Point3D(float x, float y, float z)
        {
            this.x = x; 
            this.y = y; 
            this.z = z;
        }

        public Point3D(float[] values)
        {
            if (values.Length != 3)
            {
                return;
            }
            this.x = values[0];
            this.y = values[1];
            this.z = values[2];
        }

        // Function to help inter-point comparison
        public bool Equals(Point3D otherPoint)
        {
            if(otherPoint == null)
            {
                return false;
            }
            // if comparing the same object, then they are surely equal
            if(Object.ReferenceEquals(this, otherPoint))
            {
                return true;
            }
            // only makes sense to compare apples with apples
            if(this.GetType() != otherPoint.GetType())
            {
                return false;
            }
            return this.x == otherPoint.x && this.y == otherPoint.y && this.z == otherPoint.z;
        }
    }
}