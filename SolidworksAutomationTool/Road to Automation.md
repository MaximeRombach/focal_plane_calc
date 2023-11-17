# A collection of notes for future developers

## Reference frames in Solidworks

According to an [example ](https://help.solidworks.com/2022/english/api/sldworksapi/Transform_Coordinates_from_Sketch_to_Model_Space_Example_VB.htm?verRedirect=1)buried in Solidworks' API docs, if you select a sketch point from a 3D sketch, the sketch point is already in the global (model) frame. But if the sketch point is selected from a 2D sketch, it's coordinates are in the local sketch frame.

This is a very important fact to note, as manipulating points represented in different reference frames is one of the worst things to do.

## How to know if an operation failed?

When talking directly to Solidworks using the APIs, we often have limited ways to check the status of an operation. Many of the APIs will just produce NULL if the operation failed. You won't get an error message unfortunately.




