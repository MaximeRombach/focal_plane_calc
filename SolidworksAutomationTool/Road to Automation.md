# A collection of notes for future developers

## Reference frames in Solidworks

According to an [example ](https://help.solidworks.com/2022/english/api/sldworksapi/Transform_Coordinates_from_Sketch_to_Model_Space_Example_VB.htm?verRedirect=1)buried in Solidworks' API docs, if you select a sketch point from a 3D sketch, the sketch point is already in the global (model) frame. But if the sketch point is selected from a 2D sketch, it's coordinates are in the local sketch frame.

This is a very important fact to note, as manipulating points represented in different reference frames is one of the worst things to do.

## How to know if an operation failed?

When talking directly to Solidworks using the APIs, we often have limited ways to check the status of an operation. Many of the APIs will just produce NULL if the operation failed. You won't get an error message unfortunately.

I suggest doing a NULL check before you use some variable returned from the APIs.

## Return values from the APIs

You will realize a mismatch between the data type returned by the API's and the data type specified in the Solidworks API documentation. For example, the documentation might say this API returns an array of doubles but when you hover your mouse on the API function, the return is just an object. 

This is a very common problem. I don't know why Solidworks decided to return objects in many cases. This behavior is confusing and inconsistent. 

I suggest explicitly casting the return values from an API function to the data type specified in the documentation.  

## Copying sketches will not retain the equations used in dimensions

Another problem we encountered. Solidworks does not copy equations. The current solution is to use derived sketch to have all dimensions referred to the original sketch. And since the dimensions in the original sketch are defined with global variables (equations), changes on the global variables will reflect on all sketches.

## The automation tool behaves differently from 5 mins ago

This is likely due to the miscommunication between solidworks and the automation tool. I suggest closing all solidworks instances in the Task Manager. Then relaunch the automation tool.

If you are developing the automation tool, then it's common to encounter this issue when you did a NULL dereference or other things that crashed the automation tool. The solution is the same: closing all solidworks instances in the Task Manager. 


