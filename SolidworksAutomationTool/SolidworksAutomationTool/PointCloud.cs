using System.Globalization;

namespace SolidworksAutomationTool
{
    /// <summary>
    /// This enum defines the allowed units in the txt file to be parsed in. 
    /// Developers can populate this enum to include other units in the future.
    /// </summary>
    public enum Units: ushort
    {
        Millimeter = 0,
        Meter = 1,
    }

    /// <summary>
    /// Point Cloud class that supports reading from a txt file and generating a list of 3D points out of it.
    /// </summary>
    public class PointCloud
    {
        /// <summary>
        /// A list to store all the parsed 3D points
        /// </summary>
        public List<Point3D> point3Ds = new();
        
        /// <summary>
        /// A list to store the orientaion flags of the modules whose center is colinear with the 3D points
        /// A "true" at means the corresponding module/triangle is upright. A "false" means the corresponding module/triangle is upside-down.
        /// </summary>
        public List<bool> moduleOrientations = new();

        /* The backbone of reading a point cloud from a txt file. 
         * Return:  false if the provided txt file does NOT exist / error occurred while reading the txt file.
         *          true if the txt file is read successfully. The resulting list of 3D points is stored in the class attribute point3Ds
         */

        /// <summary>
        /// Simple points-printing function. Could be used to verify the points in the text file are parsed correctly.
        /// </summary>
        public void PrintPoint3Ds()
        {
            // this number controls the field width (the fixed space taken in the console) for each number.
            const int FieldWidthRightAligned = 15;
            int pointIndex = -1;
            foreach (Point3D point in point3Ds)
            {
                pointIndex += 1;
                Console.WriteLine($"Point {pointIndex, 5}:    x: {point.x,FieldWidthRightAligned:F3} m,    y: {point.y,FieldWidthRightAligned:F3} m,    z: {point.z,FieldWidthRightAligned:F3} m");
            }
        }

        /// <summary>
        /// Simple orientation-printing function. Could be used to verify the orientation flags are parsed correctly.
        /// </summary>
        public void PrintModuleOrientations()
        {
            int moduleOrientationIndex = -1;
            foreach (bool moduleOrientation in moduleOrientations)
            {
                moduleOrientationIndex += 1;
                string orientation = moduleOrientation switch
                {
                    true => "Up",
                    false => "Down",
                };
                Console.WriteLine($"Point {moduleOrientationIndex,5}:    module orientated: {orientation}");
            }
        }

        /// <summary>
        /// Reads a txt file and parse all points into a list of 3D points
        /// </summary>
        /// <param name="fileName">The complete text file path as a string</param>
        /// <param name="dataUnit">One of the data units defined in the Units enum</param>
        /// <returns>
        /// true if the parse was sucessful.
        /// false otherwise.
        /// </returns>
        /// <exception cref="FileNotFoundException">When the path of the text file is not valid or cannot find the text file</exception>
        /// <exception cref="InvalidDataException"></exception>
        public bool ReadPointCloudFromTxt(string fileName, Units dataUnit)
        {
            // security check
            if ( !File.Exists(fileName) )
            {
                throw new FileNotFoundException($"ERROR: Cannot read point cloud from {fileName}. File doesn't exist");
            }

            if ( !Enum.IsDefined( typeof( Units ), dataUnit ) )
            {
                Console.WriteLine("ERROR: the unit entered in <function>ReadPointCloudFromTxt is not valid");
                return false;
            }

            using (StreamReader streamReader = new(fileName))
            {
                // create culture info. It will be used to specify number convention during string-to-number conversion
                CultureInfo usNumberFormatConvention = CultureInfo.CreateSpecificCulture("en-US");

                // line number starts from 0. 
                uint currentLine = 0;
                // define delimiters
                char[] splitOptions = { ' ' };

                // do a first read to go over the header. 
                string? lineRead = streamReader.ReadLine();

                // start actually reading the first line of data
                lineRead = streamReader.ReadLine() ;
                // check the number of columns in the the point cloud file
                string[] splittedLine = lineRead.Split(splitOptions, StringSplitOptions.RemoveEmptyEntries);
                // if a line contains anything other than 4 numbers, the data format is wrong. We need 3 points to define a point in 3D. 
                if (splittedLine.Length == 3)
                {
                    Console.WriteLine($"Assuming the point cloud file {fileName} has no module orientation info column");
                }
                else if (splittedLine.Length == 4)
                {
                    Console.WriteLine($"Assuming the point cloud file {fileName} has module orientation info column");
                }
                else
                {
                    throw new InvalidDataException($"ERROR: Wrong data format in line {currentLine}. This line has {splittedLine.Length} numbers instead of 3 or 4");
                }

                currentLine += 1;

                while (lineRead != null)
                {
                    // split the line read by space
                    splittedLine = lineRead.Split( splitOptions, StringSplitOptions.RemoveEmptyEntries );

                    // if a line contains anything other than 4 numbers, the data format is wrong. We need 3 points to define a point in 3D. 
                    if (!(splittedLine.Length == 3 || splittedLine.Length == 4))
                    {
                        throw new InvalidDataException($"ERROR: Wrong data format in line {currentLine}. This line has {splittedLine.Length} numbers instead of 3 or 4");
                    }

                    // Try to convert the splitted line into 3 numbers and create a Point3D instance
                    double[] convertedNumbers = new double[3];
                    for (uint axis = 0; axis < 3; axis ++)
                    {   
                        // convert strings to numbers using the US number format, a.k.a using dots as decimal points
                        if (!double.TryParse(splittedLine[axis], usNumberFormatConvention, out convertedNumbers[axis]))
                        {
                            throw new InvalidDataException($"ERROR: Wrong data format in line {currentLine} in file {fileName}. Please check for typos");
                        }

                        // Since solidworks's api seems to only take data in meters. We need to convert everything to meters, when the point cloud txt file contains data in mm
                        convertedNumbers[axis] = dataUnit switch
                        {
                            Units.Millimeter => convertedNumbers[axis] / 1000.0,
                            Units.Meter => convertedNumbers[axis],
                            // Currently only support meters and milimeters. Other units might be added later
                            _ => -1,
                        };
                    }
                    // add the create point to the Point3D list
                    point3Ds.Add( new Point3D(convertedNumbers) );

                    // store the module orientation for this "pair of points", if there are 4 columns in the point cloud text
                    if (splittedLine.Length == 4)
                    {
                        int orientationFlag = -1;
                        if (!int.TryParse(splittedLine[^1], usNumberFormatConvention, out orientationFlag) )
                        {
                            throw new InvalidDataException($"ERROR: Wrong orientation flag in line {currentLine} in file {fileName}. Please check for typos");
                        }
                        moduleOrientations.Add(orientationFlag switch
                        {
                            1 => true,  // upright triangle
                            0 => false, // upside-down triangle
                            _ => throw new InvalidDataException($"ERROR: Invalid orientation flag in line {currentLine} in file {fileName}")
                        });
                    }

                    // read the next line
                    lineRead = streamReader.ReadLine();
                    currentLine += 1;
                }
            }
            // in the case where the txt file is empty, let the user know.
            if (point3Ds.Count == 0)
            {
                throw new InvalidDataException($"WARNING: no point cloud found in {fileName}");
            }
            Console.WriteLine($"Successfully read {point3Ds.Count} points from txt");
            return true;
        }
    }

    /// <summary>
    /// A simple class to hold a 3D point.
    /// The class has 3 public double values representing the x, y and z component of a 3D point
    /// </summary>
    public class Point3D : IEquatable<Point3D>
    {
        public double x = 0.0f;
        public double y = 0.0f;
        public double z = 0.0f;

        /// <summary>
        /// Default constructors of class Point3D
        /// </summary>
        public Point3D() { }

        /// <summary>
        /// Constructor that takes 3 doubles to create a 3D point instance
        /// </summary>
        /// <param name="x"></param>
        /// <param name="y"></param>
        /// <param name="z"></param>
        public Point3D(double x, double y, double z)
        {
            this.x = x;
            this.y = y;
            this.z = z;
        }

        /// <summary>
        /// Constructor that takes an double array of 3 elements to create an instance of Point3D class.
        /// </summary>
        /// <param name="values"></param>
        /// <exception cref="InvalidDataException">Throws when the input double array does not have 3 doubles</exception>
        public Point3D(double[] values)
        {
            if (values.Length != 3)
            {
                throw new InvalidDataException("The double array used to construct a Point3D instance does not have 3 elements");
            }
            x = values[0];
            y = values[1];
            z = values[2];
        }

        /// <summary>
        /// Function to help inter-point comparison
        /// </summary>
        /// <param name="otherPoint">Instance of Point3D class</param>
        /// <returns>
        /// true if the other Point3D instance is value-wise equal to this Point3D instance
        /// false otherwise
        /// </returns>
        public bool Equals(Point3D otherPoint)
        {
            if (otherPoint == null)
            {
                return false;
            }
            // if comparing the same object, then they are surely equal
            if (Object.ReferenceEquals(this, otherPoint))
            {
                return true;
            }
            // only makes sense to compare apples with apples
            if (GetType() != otherPoint.GetType())
            {
                return false;
            }
            return x == otherPoint.x && y == otherPoint.y && z == otherPoint.z;
        }

        /// <summary>
        /// Create a mid point in between two points. The two points are passed in as "const reference"
        /// </summary>
        /// <param name="firstPoint">Constant Point3D instance</param>
        /// <param name="secondPoint">Constant Point3D instance</param>
        /// <returns>
        /// A newly created Point3D instance that spatially lives midway between the two input Point3D instances
        /// </returns>
        public static Point3D CreateMidPoint(in Point3D firstPoint, in Point3D secondPoint)
        {
            return new Point3D( firstPoint.x / 2 + secondPoint.x / 2, 
                                firstPoint.y / 2 + secondPoint.y / 2, 
                                firstPoint.z / 2 + secondPoint.z / 2 );
        }
    }

}