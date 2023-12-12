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

## Wrap repetitively used Solidworks APIs in functions

Solidworks often deprecates APIs. In order to be future proof, we try to have fewer direct API calls in the automation tool. Commonly used API calls are often wrapped in functions. 

In the case of updating the deprecated API calls, the developer will change the underlying API calls instead of the wrapper function used in the automation tool's main program.

An example of Solidworks' chain of deprecated API calls:

![SolidworksDeprecatedAPIs](https://github.com/MaximeRombach/focal_plane_calc/assets/85515041/b01dc626-42f5-4c9e-9058-afe5a538ab80)

## An API is not doing what it supposed to do

Sometimes you might realize some lastest API, let's call it "someFunction3", is not doing what it should do. 

You digged into the API documentation and found no example. However, there is a deprecated API named "someFunction2". 

Well, I suggest to give the deprecated API a try. It might do what you expected. So far I had this issue once and had to fall back to the deprecated API.

## The API does not have examples

Sometimes Solidworks only provides VBA (Visual Basic for Applications) examples or even no examples at all for some APIs.

- If there are VBA examples, you can reference the logic to write C#. The APIs in VBA are quite similar to those in C#.
- If there is no example at all, you could look online. Given that using Solidworks APIs is a niche field, your chance to find good examples online is also slim. 
  I suggest to just try the API and see if it behaves as expected. 

## Use the Macro recording function in Solidworks

A simple way to quickly create a function to achieve some part automation is to do it in Solidworks GUI manually and record the process as a macro.

Solidworks will record a whole bunch of reference frame transformations but you normally can ignore those. 

Focus on how Solidworks made your part with APIs. Once you understood the logic, it wil be easier to implement the logic in C#.

Note that this method does not always work. For example, Solidworks often does area selection when you select a few sketch segments in the GUI. 
However, area selection can be very hard to use when you select sketch segments with APIs since you need to know in which area the segments are in.

Also, some operations are not recorded in the macro. For example, the operation corresponding to flipping the normal of a plane is not recorded in the macro.
