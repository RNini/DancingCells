//Renames the open file to simplify things below

    rename("Parent")

//Sets the min and max range so the image is viewable on a computer monitor
	//Also verlays LUTs to better distinguish the channels by eye

    Stack.setChannel(1);
        setMinAndMax(100, 750);
        run("Blue");

    Stack.setChannel(2);
        setMinAndMax(800, 1800);
        run("Cyan Hot");

    Stack.setChannel(3);
        setMinAndMax(80, 1300);
        run("Red");

    Stack.setChannel(4);
        setMinAndMax(80, 1800);
        run("Thallium");

    Stack.setChannel(5);
        setMinAndMax(700, 2000);
        run("Grays");

//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

//Develops a Pop1 cell stack with the nuclear channel and Pop1 Cytoplasmic channel
	//This stack is used to track Pop1 cells
    selectWindow("Parent")
    
    Stack.setChannel(1);
        run("Reduce Dimensionality...", "frames keep");
        rename("Pop1 N1");

    selectWindow("Parent")
        Stack.setChannel(4);
        run("Reduce Dimensionality...", "frames keep");
        rename("Pop1 C1")
        setAutoThreshold("Otsu stack");
    
    imageCalculator("Min create 32-bit stack", "Pop1 N1","Pop1 C1");
        rename("Pop1")
        run("Green");
    
    close("Pop1 N1")
    close("Pop1 C1")
    
//Develops a Pop2 stack with the nuclear channel and DiI channel
	//This stack is used to track Pop2 cells
    selectWindow("Parent")

    Stack.setChannel(1);
        run("Reduce Dimensionality...", "frames keep");
        rename("Pop2 N1");

    selectWindow("Parent")
        Stack.setChannel(3);
        run("Reduce Dimensionality...", "frames keep");
        rename("Pop2 C1")
        setAutoThreshold("Otsu stack");

    imageCalculator("Min create 32-bit stack", "Pop2 N1","Pop2 C1");
        rename("Pop2 Cells")
        run("Red");
        
    close("Pop2 N1")
    close("Pop2 C1")

//Combines partitioned cells in a new window for quality control/assurance
run("Merge Channels...", "c1=[Pop1 Cells] c2=[Pop2 Cells] create keep");