function detectCraneRotationConstraint(polyhedralSurfaces, origin, radiusOfCone, heightOfCone, boomClearance, boomLength, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius,
    radiusCounterweights, heightCounterWeights, superliftBuffer, pickPoint, setPoint, minimumRadiusCone){
    // import Delaunator from 'delaunator';
    import ('d3').then(d3 => {
        function distMeshCone (origin, radius, height) {
            // pfix is a fixed set of points, which are given as input to the system.
            // The geometry is given as a distance function fd. This function returns the signed distance from each node location p to the closest boundary.
            // Originally, fd and fh are considered wxternal functions. But in this case, I am making them implicit. 
            // The (relative) desired edge length function h(x, y) is given as a function fh, which returns h for all input points.
    
            // Starting with a uniform grid. I think changing the grid is a small thing. It can be easily done later on. 
            function fdCircle (points, radius) {
                var boundaryDistance = Math.sqrt(points[0]*points[0] + points[1]*points[1])-radius;
                return boundaryDistance;
            }
    
            // To start with, I am considering a constant value of h for a uniform meshing. For non uniform meshing as well, there can be specific functions generated. 
            function fhCircle (points) {
                var edgeLength = 2;
                return edgeLength;
            }
    
            var inputh0 = 2; 
            var inputRadiusCone = radius; 
            var inputConeHeight = height; 
            var inputBboxCone = [
                [-inputRadiusCone,-inputRadiusCone],
                [inputRadiusCone,inputRadiusCone]
            ];
            var inputPfixCone = [
                [0,0], 
                [inputRadiusCone,0], 
                [0,inputRadiusCone],
                [-inputRadiusCone,0], 
                [0,-inputRadiusCone]
            ];
            var dptol=.000001; 
            var ttol=.1; 
            var Fscale=1.2; 
            var deltat=.2; 
            var geps=.000001*inputh0; 
            var deps=Math.sqrt(geps)*inputh0;
    
            // Creating an initial meshgrid
            const x_min = inputBboxCone[0][0];
            const y_min = inputBboxCone[0][1];
            const x_max = inputBboxCone[1][0];
            const y_max = inputBboxCone[1][1];
    
            const rows_x = Math.floor((x_max - x_min) / inputh0) + 1;
            const rows_y = Math.floor((y_max - y_min) / (inputh0 * Math.sqrt(3) / 2)) + 1;
            // rows_x and rows_y will contain the number of rows and columns needed in the provided bounding box dimensions based on the given value of h0. 
    
            const x = [];
            const y = [];
    
            for (let i = 0; i < rows_y; i++) {
                const row_x = [];
                const row_y = [];
            
                for (let j = 0; j < rows_x; j++) {
                    row_x.push(x_min + j * inputh0);
                    row_y.push(y_min + i * inputh0 * Math.sqrt(3) / 2);
                }
            
                x.push(row_x);
                y.push(row_y);
            }
            // From here, x and y will be an array of arrays. For the corresponding addresses, each elements in the x and y will be an array. 
            // The elements of the array within x and y will point to the corresponding coordinates. 
            // For instance, the array at the [0] index in both x and y will provide the x and y values corresponding to the coordinates at the first row of points. 
            
            // Shifting the x coordinates in the alternate rows by h0/2. 
            for (let i = 1; i < x.length; i += 2) {
                for (let j = 0; j < x[i].length; j++) {
                    x[i][j] += inputh0/2;
                }
            }
    
            // The .flat() functions converts the x and y from an array of arrays into a single array. 
            // All the elements are stacked in a series. 
            const flattenedX = x.flat();
            const flattenedY = y.flat();
            // console.log(flattenedX);
            //console.log(flattenedY.length);
    
            // Combine flattenedX and flattenedY column-wise
            // p developed from here will be an array of arrays, where each inside array will be a coordinate. 
            var p = flattenedX.map((value, index) => [value, flattenedY[index]]);
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // console.log(p.length);
            var filteredP = [];
    
            // Removing the points outside the boundary.
            // Only the interior points with negative distances (allowing a tolerance geps) are kept. 
            for (let i=0; i<p.length; i++) {
                const row = p[i];
    
                if (fdCircle(row, radius) < geps) {
                    filteredP.push(row);
                }
            }
            // Again, the filteredP will also be same as P in terms of the type. It will contain the filtered points.  
            // console.log(filteredP.length);
    
            // Assigning filteredP back to p because all the further operations are done on p. 
            p = [...filteredP];
            // console.log(p.length);
            const r0 = [];
    
            // As per the paper, r0 is showing the probability to keep a point. So, it most probably will be just a value. 
            for (let i = 0; i < p.length; i++) {
                const row = p[i];
                const evaluation = fhCircle(row);
                const inverseSquare = 1 / evaluation ** 2;
                r0.push(inverseSquare);
            }
    
            // Generating a random array of the length of p.length. 
            const randValues = Array.from({ length: p.length }, () => Math.random());
            // console.log(randValues[0]);
    
            // Normalize the values in r0. 
            const normalizedR0 = r0.map(value => value / Math.max(...r0));
    
            // Selecting rows based on the condition
            const selectedRows = [];
    
            // I think, logically, this is doing nothing. I can probably remove this as well. What difference does the comparsion of a random value with its normalized value make? I am not sure with this. 
            for (let i = 0; i < p.length; i++) {
                if (randValues[i] < normalizedR0[i]) {
                selectedRows.push(p[i]);
                }
            }
    
            // Concatenating matrices. So, what I understand is, concatenatedP will also be an array of arrays, where each point coordinates are represented.
            const concatenatedP = inputPfixCone.concat(selectedRows);
            p = [...concatenatedP];
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // Storing the length of p array in the variable N. 
            const N = p.length;
            //console.log(N);
    
            // Initializing pold with infinity
            let pold = Array.from({ length: p.length }, () => [Infinity, Infinity]);
            var numIter = 1;
            while (true) {
                // 3. Retriangulation by the Delaunay algorithm
                var displacement = 0;
                var difference = [];
                //Computing the elementwise difference between the p and pold. 
                for (var k=0; k<p.length; k++) {
                    difference[k] = [];
                    difference[k][0] = p[k][0] - pold[k][0];
                    difference[k][1] = p[k][1] - pold[k][1];
                }
    
                // Computing the square root of sum of squares of the difference and dividing by h0.
                var normalizedDifference = []; 
                for (var l=0; l<difference.length; l++) {
                    normalizedDifference[l] = Math.sqrt(difference[l][0]*difference[l][0] + difference[l][1]*difference[l][1])/inputh0;
                }
                // Finding the maximum value of the normalized difference from the array. 
                displacement = Math.max(...normalizedDifference);
                // console.log(displacement);
                //p.map((point, index) => Math.sqrt(point.reduce((sum, coord, axis) => sum + Math.pow(coord - pold[index][axis], 2), 0)) / h0);
            
                if (displacement > ttol) {
                    pold = p.map(point => [...point]); // Save current positions (deep copy). 
                    // We need to do this because if we directly assign pold = p, both will reference to the same object. So making changes in one would afftect the other as well. 
                    // console.log(pold[0]);
                    // So, since I am using the Delaunator library from the Javascript, the results are accordingly. 
                    // Delaunator.from(p).triangles will return an array. 
                    // Each consecutive three values in the triangles will show the indices of the points in p, which form a triangle in Delaunay.
                    // const trianglesIndices = Delaunator.from(p).triangles;
                    const trianglesIndices = new d3.Delaunay(p.flat()).triangles;
                    // console.log(trianglesIndices);
    
                    // Create an array of triangle points similar to what was needed in the original code, so that it can be processed accordingly.
                    const triangles = Array.from({ length: trianglesIndices.length/3 }, () => Array(3));
                    for (let i=0; i<trianglesIndices.length; i+=3){
                        for (let j=0; j<3; j++){
                            triangles[i/3][j] = trianglesIndices[i+j];
                        }
                    }
                    // console.log(triangles.length);
    
    
                    // Compute centroids of all the created triangles
                    const pmid = triangles.map(t => [
                        (p[t[0]][0] + p[t[1]][0] + p[t[2]][0]) / 3,
                        (p[t[0]][1] + p[t[1]][1] + p[t[2]][1]) / 3
                    ]);
                    // console.log(pmid[0]);
            
                    // Keep interior triangles. It filters those triangles whose centeroids are inside the boundary region. 
                    // Filter method will keep the triangles for which the condition is true. 
                    var interiorTriangles = triangles.filter((t, index) => fdCircle(pmid[index], radius) < -geps); //Although the paper uses -geps here, I think it should be geps. 
                    // Looking into more detail, it does not make much difference geps or -geps). 
                    // console.log(interiorTriangles.length);
            
                    // Describe each bar by a unique pair of nodes. So basically from the selected nodes, you are creating the bars. 
                    var bars = interiorTriangles.flatMap(t => [
                        [t[0], t[1]],
                        [t[0], t[2]],
                        [t[1], t[2]]
                    ]);
                    // console.log(bars);
                    // To obtain a unique set of sorted bars from the bars saved in the previous step. 
                    const sortedBars = bars.map(bar => [...bar].sort((a, b) => a - b));
                    const uniqueBars = Array.from(new Set(sortedBars.map(JSON.stringify)), JSON.parse);
                    // Assigning uniqueBars back to the original bars, which is the variable used ahead in the code. 
                    bars = [...uniqueBars];
                    // console.log(bars);
                    // Graphical output of the current mesh (console.log used for demonstration)
                    //console.log("Iteration:", iterationCount);
                    //console.log("Triangles:", interiorTriangles);
                    //console.log("Bars:", uniqueBars);
            
                    //iterationCount++;
            
                    // Additional code for graphical output (adjust as needed)
                }
                // Now on, you will try to move the bar points based on the bar lengths and forces. 
                // So now, the indices stored in the bars need to be taken into the form of vector representing the difference between x and y components. 
                // List of bar vectors
                const barvec = bars.map(bar => [
                    p[bar[0]][0] - p[bar[1]][0],
                    p[bar[0]][1] - p[bar[1]][1]
                    ]);
                // console.log(barvec);
                // To my understanding, barvec will represent the x and y components of the vector representing the bars. 
    
    
                // Compute the length of each vector in the barvec.
                const L = barvec.map(vector => Math.sqrt(vector[0] ** 2 + vector[1] ** 2));
                // console.log(L);
    
                // Computing the midpoint of all the edges and then computing the fh function for these mid-point coordinates. 
                // So, hbars will be an array containing the function values of fh evaluated at the midpoints of the bars. 
                // Each element of hbars corresponds to a different bar.
                // Compute hbars
                const hbars = bars.map(bar => {
                    const midpoint = [
                    (p[bar[0]][0] + p[bar[1]][0]) / 2,
                    (p[bar[0]][1] + p[bar[1]][1]) / 2
                    ];
                    return fhCircle(midpoint);
                });
    
                // Computing the desired length of the bars
                const sumL2 = L.reduce((sum, value) => sum + value ** 2, 0);
                const sumHbars2 = hbars.reduce((sum, value) => sum + value ** 2, 0);
    
                const L0 = hbars.map((hbar) => hbar * Fscale * Math.sqrt(sumL2 / sumHbars2));
                // console.log(L0);
                // Computing the bar forces. It is being done by computing the difference between desired length and the computed length. 
                const F = L0.map((L0i, index) => Math.max(L0i - L[index], 0));
                // console.log (F);
                // Bar forces (x, y components)
                const Fvec = F.map((force, index) => [
                    force / L[index] * barvec[index][0],
                    force / L[index] * barvec[index][1]
                ]);
                // console.log(Fvec);
                //The resulting Fvec array contains arrays where each element represents the force vector (x and y components) for each bar.
                // console.log(Fvec);
                // Create a matrix representing the summation of forces at each node. 
    
                // var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                // The length of Ftot is same as that of the length of p (storing total number of nodes). This is to my understanding. 
    
                // Ftot has to be the sum of force vectors in Fvec, from all the bars meeting at a node.
                // A stretching force has a positive sign, and its direction is given by the two component vectors in bars. I think there is something wrong in this.  
                /*for (var z=0; z < p.length; z++) {
                    for (var a = 0; a < bars.length; a++) {
                        if (bars[a].includes(z)) {
                            Ftot[z][0] = Ftot[z][0] + Fvec[a][0];
                            Ftot[z][1] = Ftot[z][1] + Fvec[a][1];
                        }
                    }
                }*/
                // console.log(Ftot[0].length);
    
                var sparseInput1 = [];
                for (var i = 0; i < bars.length; i++) {
                    sparseInput1[i] = [];
                    sparseInput1[i][0] = bars[i][0];
                    sparseInput1[i][1] = bars[i][0];
                    sparseInput1[i][2] = bars[i][1];
                    sparseInput1[i][3] = bars[i][1];
                }
                var flatSparseInput1 = sparseInput1.flat();
    
                var originalArray = [0, 1, 0, 1];
                var lenF = F.length;
                var sparseInput2 = Array.from({ length: F.length }, () => [...originalArray]);
                var flatSparseInput2 = sparseInput2.flat();
    
                var sparseInput3 = [];
    
                for (var i = 0; i < Fvec.length; i++) {
                    sparseInput3[i] = [];
                    sparseInput3[i][0] = Fvec[i][0];
                    sparseInput3[i][1] = Fvec[i][1];
                    sparseInput3[i][2] = -1*Fvec[i][0];
                    sparseInput3[i][3] = -1*Fvec[i][1];
                }
                var flatSparseInput3 = sparseInput3.flat();
    
                var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                for (var i = 0; i < flatSparseInput1.length; i++) {
                    Ftot[flatSparseInput1[i]][flatSparseInput2[i]] += flatSparseInput3[i];
                }
    
    
                //bars.forEach((bar, index) => {
                //    Ftot[bar[0]][0] += Fvec[index][0];
                //    Ftot[bar[0]][1] += Fvec[index][1];
                //    Ftot[bar[1]][0] -= Fvec[index][0];
                //    Ftot[bar[1]][1] -= Fvec[index][1];
                //});
    
                // Setting up force components to be zero at the fixed nodes. 
                // Force = 0 at fixed points
                const numRows = inputPfixCone.length; // Assuming pfix is a 2D array
                const numCols = Ftot[0].length; // Assuming Ftot is a 2D array
                for (let i = 0; i < numRows; i++) {
                    for (let j = 0; j < numCols; j++) {
                    Ftot[i][j] = 0;
                    }
                }
                //pfix.forEach(fixedPointIndex => {
                //    Ftot[fixedPointIndex][0] = 0;
                //    Ftot[fixedPointIndex][1] = 0;
                //});
    
                // Update node positions
                p = p.map((point, index) => [
                    point[0] + deltat * Ftot[index][0],
                    point[1] + deltat * Ftot[index][1]
                ]);
    
                // Finds points in p which are outside the d>0.
                // Resulting ix will be a boolean array, containing true if the distance > 0, false otherwise. 
                const d = p.map(point => fdCircle(point, radius));
                const ix = d.map(value => value > 0);
    
                // For the points where the value of d>0, you want to move them inwards. 
                // The first step is to find out the gradient of those points.
                // Then, based on the gradient, the points are moved. 
                const projectedPoints = ix.map((value, index) => {
                    if (value) {
                        // Here you make slight changes in the coordinate values at p. 
                        const perturbedX = [p[index][0] + deps, p[index][1]];
                        const perturbedY = [p[index][0], p[index][1] + deps];
                
                        // Here, you are trying to find out the gradient I guess. To see how much change should you do in the x and y coordinate. 
                        const dgradx = (fdCircle(perturbedX, radius) - d[index]) / deps;
                        const dgrady = (fdCircle(perturbedY, radius) - d[index]) / deps;
                        
                        // Here, you are making those changes in the x and y coordinate values. 
                        return [
                            p[index][0] - d[index] * dgradx,
                            p[index][1] - d[index] * dgrady
                        ];
                    } else {
                        return p[index];
                    }
                });
                // So the projected points here will be an array containing all the points in p, which also include the points which have been moved inwards. 
                // Update p with the projected points
                for (let i = 0; i < p.length; i++) {
                    p[i] = projectedPoints[i];
                }
    
                // Breaking criteria to come out of the loop. 
                // It checks if all the interior nodes move less than dptol. 
                // This is running infinitely long, I don't understand why. I am putting another break condition on the maximum number of iterations as 1000 based on https://github.com/ionhandshaker/distmesh/blob/master/README.md.
                const threshold = dptol;
                const logicalArray2 = d.map(element => element < -geps); // Create a logical array based on d for which the values are less than -geps. 
                var ftotLessThanGeps = [];
                for (let i = 0; i < logicalArray2.length; i++) {
                    if (logicalArray2[i]) {
                        ftotLessThanGeps.push(Ftot[i]);
                    }
                }
                var squaredNorms = [];
                for (let i = 0; i < ftotLessThanGeps.length; i++) {
                    squaredNorms[i] = Math.sqrt((ftotLessThanGeps[i][0]*deltat)**2 + (ftotLessThanGeps[i][1]*deltat)**2);
                }
                const maxSquaredNorm = Math.max(...squaredNorms);
                // const maxNorm = Math.sqrt(maxSquaredNorm);
                // console.log(maxSquaredNorm / h0);
                if (maxSquaredNorm / inputh0 < threshold) {
                    // return {p, bars}; 
                    break;
                } /*else if (numIter === 2000) {
                    return {p, bars}; 
                    break;
                } else {
                    numIter = numIter+1;
                }*/
                // Additional conditions and iterations may be added here
            }
            var tanTheta = height/radius ;
            var conePoints = [];
            for (var i=0; i<p.length; i++) {
            conePoints[i] = [];
            conePoints[i][0] = p[i][0];
            conePoints[i][1] = p[i][1];
            conePoints[i][2] = Math.sqrt(p[i][0]*p[i][0]+p[i][1]*p[i][1])*tanTheta;
            }
            // Transforming the coordinates. 
            for (var i = 0;  i < conePoints.length; i++) { 
                conePoints[i][0] += origin[0];
                conePoints[i][1] += origin[1];
                conePoints[i][2] += origin[2];
            }
            var finalPoints = [];
            for (let i=0; i<conePoints.length; i++) {
                finalPoints[i] = [];
                finalPoints[i][0] = parseFloat(conePoints[i][0].toFixed(6));
                finalPoints[i][1] = parseFloat(conePoints[i][1].toFixed(6));
                finalPoints[i][2] = parseFloat(conePoints[i][2].toFixed(6));
            }
            return {finalPoints, bars, interiorTriangles};
        }
    
        //var radiusCircle = 10;
        //var heightCone = 5; 
        //var origin = [0,0,0]; 
        //var triangulatedCone = distMeshCone (origin, radiusCircle, heightCone);
        // console.log(triangulatedCircle.bars)
        /*for (var m = 0; m<triangulatedCone.p.length; m++) {
            const subArray = triangulatedCone.p[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCone.bars.length; m++) {
            const subArray = triangulatedCone.bars[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCone.finalPoints.length; m++) {
            const subArray = triangulatedCone.finalPoints[m];
            console.log(subArray);
        }*/
    
        function distMeshCylinder (origin, radius, height) {
            // pfix is a fixed set of points, which are given as input to the system.
            // The geometry is given as a distance function fd. This function returns the signed distance from each node location p to the closest boundary.
            // The (relative) desired edge length function h(x, y) is given as a function fh, which returns h for all input points.
            function fdRectangle (points, boundingBox) {
                var xdistance = Math.min(-1*boundingBox[0][0] + points[0], boundingBox[1][0] - points[0]);
                var ydistance = Math.min(-1*boundingBox[0][1] + points[1], boundingBox[1][1] - points[1]);
                var boundaryDistance = -1 * Math.min(xdistance, ydistance);
                return boundaryDistance;
            }
            function fhCylinder (points) {
                var edgeLength = 1;
                return edgeLength;
            }
    
            var inputh0 = 1;
            
            var inputRadiusCylinder = radius;
    
            var inputCylinderHeight = height; 
    
            var inputPfixCylinder = [
                [-Math.PI*inputRadiusCylinder,-inputCylinderHeight/2], 
                [Math.PI*inputRadiusCylinder,-inputCylinderHeight/2], 
                [Math.PI*inputRadiusCylinder,inputCylinderHeight/2],
                [-Math.PI*inputRadiusCylinder,inputCylinderHeight/2]
            ];
            
            var inputBboxCylinder = [
                [-Math.PI*inputRadiusCylinder,-inputCylinderHeight/2],
                [Math.PI*inputRadiusCylinder,inputCylinderHeight/2]
            ];
    
            var dptol=.001; 
            var ttol=.1; 
            var Fscale=1.2; 
            var deltat=.2; 
            var geps=.001*inputh0; 
            var deps=Math.sqrt(geps)*inputh0;
    
            // Creating an initial meshgrid
            const x_min = inputBboxCylinder[0][0];
            const y_min = inputBboxCylinder[0][1];
            const x_max = inputBboxCylinder[1][0];
            const y_max = inputBboxCylinder[1][1];
    
            const rows_x = Math.floor((x_max - x_min) / inputh0) + 1;
            const rows_y = Math.floor((y_max - y_min) / (inputh0 * Math.sqrt(3) / 2)) + 1;
            // rows_x and rows_y will contain the number of rows and columns needed in the provided bounding box dimensions based on the given value of h0. 
    
            const x = [];
            const y = [];
    
            for (let i = 0; i < rows_y; i++) {
                const row_x = [];
                const row_y = [];
            
                for (let j = 0; j < rows_x; j++) {
                    row_x.push(x_min + j * inputh0);
                    row_y.push(y_min + i * inputh0 * Math.sqrt(3) / 2);
                }
            
                x.push(row_x);
                y.push(row_y);
            }
            // From here, x and y will be an array of arrays. For the corresponding addresses, each elements in the x and y will be an array. 
            // The elements of the array within x and y will point to the corresponding coordinates. 
            // For instance, the array at the [0] index in both x and y will provide the x and y values corresponding to the coordinates at the first row of points. 
            
            // Shifting the x coordinates in the alternate rows by h0/2. 
            for (let i = 1; i < x.length; i += 2) {
                for (let j = 0; j < x[i].length; j++) {
                    x[i][j] += inputh0/2;
                }
            }
    
            // The .flat() functions converts the x and y from an array of arrays into a single array. 
            // All the elements are stacked in a series. 
            const flattenedX = x.flat();
            const flattenedY = y.flat();
            // console.log(flattenedX);
            //console.log(flattenedY.length);
    
            // Combine flattenedX and flattenedY column-wise
            // p developed from here will be an array of arrays, where each inside array will be a coordinate. 
            var p = flattenedX.map((value, index) => [value, flattenedY[index]]);
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // console.log(p.length);
            var filteredP = [];
    
            // Removing the points outside the boundary.
            // Only the interior points with negative distances (allowing a tolerance geps) are kept. 
            for (let i=0; i<p.length; i++) {
                const row = p[i];
    
                if (fdRectangle(row, inputBboxCylinder) < geps) {
                    filteredP.push(row);
                }
            }
            // Again, the filteredP will also be same as P in terms of the type. It will contain the filtered points.  
            // console.log(filteredP.length);
    
            // Assigning filteredP back to p because all the further operations are done on p. 
            p = [...filteredP];
            // console.log(p.length);
            const r0 = [];
    
            // As per the paper, r0 is showing the probability to keep a point. So, it most probably will be just a value. 
            for (let i = 0; i < p.length; i++) {
                const row = p[i];
                const evaluation = fhCylinder(row);
                const inverseSquare = 1 / evaluation ** 2;
                r0.push(inverseSquare);
            }
    
            // Generating a random array of the length of p.length. 
            const randValues = Array.from({ length: p.length }, () => Math.random());
            // console.log(randValues[0]);
    
            // Normalize the values in r0. 
            const normalizedR0 = r0.map(value => value / Math.max(...r0));
    
            // Selecting rows based on the condition
            const selectedRows = [];
    
            // I think, logically, this is doing nothing. I can probably remove this as well. What difference does the comparsion of a random value with its normalized value make? I am not sure with this. 
            for (let i = 0; i < p.length; i++) {
                if (randValues[i] < normalizedR0[i]) {
                selectedRows.push(p[i]);
                }
            }
    
            // Concatenating matrices. So, what I understand is, concatenatedP will also be an array of arrays, where each point coordinates are represented.
            const concatenatedP = inputPfixCylinder.concat(selectedRows);
            p = [...concatenatedP];
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // Storing the length of p array in the variable N. 
            const N = p.length;
            //console.log(N);
    
            // Initializing pold with infinity
            let pold = Array.from({ length: p.length }, () => [Infinity, Infinity]);
            var numIter = 1;
            while (true) {
                // 3. Retriangulation by the Delaunay algorithm
                var displacement = 0;
                var difference = [];
                //Computing the elementwise difference between the p and pold. 
                for (var k=0; k<p.length; k++) {
                    difference[k] = [];
                    difference[k][0] = p[k][0] - pold[k][0];
                    difference[k][1] = p[k][1] - pold[k][1];
                }
    
                // Computing the square root of sum of squares of the difference and dividing by h0.
                var normalizedDifference = []; 
                for (var l=0; l<difference.length; l++) {
                    normalizedDifference[l] = Math.sqrt(difference[l][0]*difference[l][0] + difference[l][1]*difference[l][1])/inputh0;
                }
                // Finding the maximum value of the normalized difference from the array. 
                displacement = Math.max(...normalizedDifference);
                // console.log(displacement);
                //p.map((point, index) => Math.sqrt(point.reduce((sum, coord, axis) => sum + Math.pow(coord - pold[index][axis], 2), 0)) / h0);
            
                if (displacement > ttol) {
                    pold = p.map(point => [...point]); // Save current positions (deep copy). 
                    // We need to do this because if we directly assign pold = p, both will reference to the same object. So making changes in one would afftect the other as well. 
                    // console.log(pold[0]);
                    // So, since I am using the Delaunator library from the Javascript, the results are accordingly. 
                    // Delaunator.from(p).triangles will return an array. 
                    // Each consecutive three values in the triangles will show the indices of the points in p, which form a triangle in Delaunay.
                    // const trianglesIndices = Delaunator.from(p).triangles;
                    const trianglesIndices = new d3.Delaunay(p.flat()).triangles;
                    // console.log(trianglesIndices);
    
                    // Create an array of triangle points similar to what was needed in the original code, so that it can be processed accordingly.
                    const triangles = Array.from({ length: trianglesIndices.length/3 }, () => Array(3));
                    for (let i=0; i<trianglesIndices.length; i+=3){
                        for (let j=0; j<3; j++){
                            triangles[i/3][j] = trianglesIndices[i+j];
                        }
                    }
                    // console.log(triangles.length);
    
    
                    // Compute centroids of all the created triangles
                    const pmid = triangles.map(t => [
                        (p[t[0]][0] + p[t[1]][0] + p[t[2]][0]) / 3,
                        (p[t[0]][1] + p[t[1]][1] + p[t[2]][1]) / 3
                    ]);
                    // console.log(pmid[0]);
            
                    // Keep interior triangles. It filters those triangles whose centeroids are inside the boundary region. 
                    // Filter method will keep the triangles for which the condition is true. 
                    var interiorTriangles = triangles.filter((t, index) => fdRectangle(pmid[index], inputBboxCylinder) < -geps); //Although the paper uses -geps here, I think it should be geps. 
                    // Looking into more detail, it does not make much difference geps or -geps). 
                    // console.log(interiorTriangles.length);
            
                    // Describe each bar by a unique pair of nodes. So basically from the selected nodes, you are creating the bars. 
                    var bars = interiorTriangles.flatMap(t => [
                        [t[0], t[1]],
                        [t[0], t[2]],
                        [t[1], t[2]]
                    ]);
                    // console.log(bars);
                    // To obtain a unique set of sorted bars from the bars saved in the previous step. 
                    const sortedBars = bars.map(bar => [...bar].sort((a, b) => a - b));
                    const uniqueBars = Array.from(new Set(sortedBars.map(JSON.stringify)), JSON.parse);
                    // Assigning uniqueBars back to the original bars, which is the variable used ahead in the code. 
                    bars = [...uniqueBars];
                    // console.log(bars);
                    // Graphical output of the current mesh (console.log used for demonstration)
                    //console.log("Iteration:", iterationCount);
                    //console.log("Triangles:", interiorTriangles);
                    //console.log("Bars:", uniqueBars);
            
                    //iterationCount++;
            
                    // Additional code for graphical output (adjust as needed)
                }
                // Now on, you will try to move the bar points based on the bar lengths and forces. 
                // So now, the indices stored in the bars need to be taken into the form of vector representing the difference between x and y components. 
                // List of bar vectors
                const barvec = bars.map(bar => [
                    p[bar[0]][0] - p[bar[1]][0],
                    p[bar[0]][1] - p[bar[1]][1]
                    ]);
                // console.log(barvec);
                // To my understanding, barvec will represent the x and y components of the vector representing the bars. 
    
    
                // Compute the length of each vector in the barvec.
                const L = barvec.map(vector => Math.sqrt(vector[0] ** 2 + vector[1] ** 2));
                // console.log(L);
    
                // Computing the midpoint of all the edges and then computing the fh function for these mid-point coordinates. 
                // So, hbars will be an array containing the function values of fh evaluated at the midpoints of the bars. 
                // Each element of hbars corresponds to a different bar.
                // Compute hbars
                const hbars = bars.map(bar => {
                    const midpoint = [
                    (p[bar[0]][0] + p[bar[1]][0]) / 2,
                    (p[bar[0]][1] + p[bar[1]][1]) / 2
                    ];
                    return fhCylinder(midpoint);
                });
    
                // Computing the desired length of the bars
                const sumL2 = L.reduce((sum, value) => sum + value ** 2, 0);
                const sumHbars2 = hbars.reduce((sum, value) => sum + value ** 2, 0);
    
                const L0 = hbars.map((hbar) => hbar * Fscale * Math.sqrt(sumL2 / sumHbars2));
                // console.log(L0);
                // Computing the bar forces. It is being done by computing the difference between desired length and the computed length. 
                const F = L0.map((L0i, index) => Math.max(L0i - L[index], 0));
                // console.log (F);
                // Bar forces (x, y components)
                const Fvec = F.map((force, index) => [
                    force / L[index] * barvec[index][0],
                    force / L[index] * barvec[index][1]
                ]);
                // console.log(Fvec);
                //The resulting Fvec array contains arrays where each element represents the force vector (x and y components) for each bar.
                // console.log(Fvec);
                // Create a matrix representing the summation of forces at each node. 
    
                // var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                // The length of Ftot is same as that of the length of p (storing total number of nodes). This is to my understanding. 
    
                // Ftot has to be the sum of force vectors in Fvec, from all the bars meeting at a node.
                // A stretching force has a positive sign, and its direction is given by the two component vectors in bars. I think there is something wrong in this.  
                /*for (var z=0; z < p.length; z++) {
                    for (var a = 0; a < bars.length; a++) {
                        if (bars[a].includes(z)) {
                            Ftot[z][0] = Ftot[z][0] + Fvec[a][0];
                            Ftot[z][1] = Ftot[z][1] + Fvec[a][1];
                        }
                    }
                }*/
                // console.log(Ftot[0].length);
    
                var sparseInput1 = [];
                for (var i = 0; i < bars.length; i++) {
                    sparseInput1[i] = [];
                    sparseInput1[i][0] = bars[i][0];
                    sparseInput1[i][1] = bars[i][0];
                    sparseInput1[i][2] = bars[i][1];
                    sparseInput1[i][3] = bars[i][1];
                }
                var flatSparseInput1 = sparseInput1.flat();
    
                var originalArray = [0, 1, 0, 1];
                var lenF = F.length;
                var sparseInput2 = Array.from({ length: F.length }, () => [...originalArray]);
                var flatSparseInput2 = sparseInput2.flat();
    
                var sparseInput3 = [];
    
                for (var i = 0; i < Fvec.length; i++) {
                    sparseInput3[i] = [];
                    sparseInput3[i][0] = Fvec[i][0];
                    sparseInput3[i][1] = Fvec[i][1];
                    sparseInput3[i][2] = -1*Fvec[i][0];
                    sparseInput3[i][3] = -1*Fvec[i][1];
                }
                var flatSparseInput3 = sparseInput3.flat();
    
                var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                for (var i = 0; i < flatSparseInput1.length; i++) {
                    Ftot[flatSparseInput1[i]][flatSparseInput2[i]] += flatSparseInput3[i];
                }
    
    
                //bars.forEach((bar, index) => {
                //    Ftot[bar[0]][0] += Fvec[index][0];
                //    Ftot[bar[0]][1] += Fvec[index][1];
                //    Ftot[bar[1]][0] -= Fvec[index][0];
                //    Ftot[bar[1]][1] -= Fvec[index][1];
                //});
    
                // Setting up force components to be zero at the fixed nodes. 
                // Force = 0 at fixed points
                const numRows = inputPfixCylinder.length; // Assuming pfix is a 2D array
                const numCols = Ftot[0].length; // Assuming Ftot is a 2D array
                for (let i = 0; i < numRows; i++) {
                    for (let j = 0; j < numCols; j++) {
                    Ftot[i][j] = 0;
                    }
                }
                //pfix.forEach(fixedPointIndex => {
                //    Ftot[fixedPointIndex][0] = 0;
                //    Ftot[fixedPointIndex][1] = 0;
                //});
    
                // Update node positions
                p = p.map((point, index) => [
                    point[0] + deltat * Ftot[index][0],
                    point[1] + deltat * Ftot[index][1]
                ]);
    
                // Finds points in p which are outside the d>0.
                // Resulting ix will be a boolean array, containing true if the distance > 0, false otherwise. 
                const d = p.map(point => fdRectangle(point, inputBboxCylinder));
                const ix = d.map(value => value > 0);
    
                // For the points where the value of d>0, you want to move them inwards. 
                // The first step is to find out the gradient of those points.
                // Then, based on the gradient, the points are moved. 
                const projectedPoints = ix.map((value, index) => {
                    if (value) {
                        // Here you make slight changes in the coordinate values at p. 
                        const perturbedX = [p[index][0] + deps, p[index][1]];
                        const perturbedY = [p[index][0], p[index][1] + deps];
                
                        // Here, you are trying to find out the gradient I guess. To see how much change should you do in the x and y coordinate. 
                        const dgradx = (fdRectangle(perturbedX, inputBboxCylinder) - d[index]) / deps;
                        const dgrady = (fdRectangle(perturbedY, inputBboxCylinder) - d[index]) / deps;
                        
                        // Here, you are making those changes in the x and y coordinate values. 
                        return [
                            p[index][0] - d[index] * dgradx,
                            p[index][1] - d[index] * dgrady
                        ];
                    } else {
                        return p[index];
                    }
                });
                // So the projected points here will be an array containing all the points in p, which also include the points which have been moved inwards. 
                // Update p with the projected points
                for (let i = 0; i < p.length; i++) {
                    p[i] = projectedPoints[i];
                }
    
                // Breaking criteria to come out of the loop. 
                // It checks if all the interior nodes move less than dptol. 
                // This is running infinitely long, I don't understand why. I am putting another break condition on the maximum number of iterations as 1000 based on https://github.com/ionhandshaker/distmesh/blob/master/README.md.
                const threshold = dptol;
                const logicalArray2 = d.map(element => element < -geps); // Create a logical array based on d for which the values are less than -geps. 
                var ftotLessThanGeps = [];
                for (let i = 0; i < logicalArray2.length; i++) {
                    if (logicalArray2[i]) {
                        ftotLessThanGeps.push(Ftot[i]);
                    }
                }
                var squaredNorms = [];
                for (let i = 0; i < ftotLessThanGeps.length; i++) {
                    squaredNorms[i] = Math.sqrt((ftotLessThanGeps[i][0]*deltat)**2 + (ftotLessThanGeps[i][1]*deltat)**2);
                }
                const maxSquaredNorm = Math.max(...squaredNorms);
                // const maxNorm = Math.sqrt(maxSquaredNorm);
                // console.log(maxSquaredNorm / h0);
                if (maxSquaredNorm / inputh0 < threshold) {
                    // return {p, bars}; 
                    break;
                } /*else if (numIter === 2000) {
                    return {p, bars}; 
                    break;
                } else {
                    numIter = numIter+1;
                }*/
                // Additional conditions and iterations may be added here
            }
            var xxx = [];
            // Getting all the X-coordinate values from the points. If the input rectangle is considered as the origin being at the centre, therefore, we transform x  to x+pie*R.  
            for (var i=0; i<p.length; i++) {
                xxx.push(p[i][0]+Math.PI*inputRadiusCylinder);
            }
            var cylinderPoints = [];
    
            for (var i=0; i<xxx.length; i++) {
            cylinderPoints[i] = [];
            cylinderPoints[i][0] = inputRadiusCylinder * Math.sin((xxx[i]/Math.max(...xxx)) * Math.PI * 2);
            cylinderPoints[i][1] = inputRadiusCylinder * Math.cos((xxx[i]/Math.max(...xxx)) * Math.PI * 2);
            cylinderPoints[i][2] = p[i][1] + inputCylinderHeight/2;
            }
            // Transforming the coordinates
            for (var i = 0;  i < cylinderPoints.length; i++) { 
                cylinderPoints[i][0] += origin[0];
                cylinderPoints[i][1] += origin[1];
                cylinderPoints[i][2] += origin[2];
            }
            var finalPoints = [];
            for (let i=0; i<cylinderPoints.length; i++) {
                finalPoints[i] = [];
                finalPoints[i][0] = parseFloat(cylinderPoints[i][0].toFixed(6));
                finalPoints[i][1] = parseFloat(cylinderPoints[i][1].toFixed(6));
                finalPoints[i][2] = parseFloat(cylinderPoints[i][2].toFixed(6));
            }
            return {finalPoints, bars, interiorTriangles};
        }
    
        //var radiusCircle = 4;
        //var heightCylinder = 3; 
        //var origin = [0,0,0]; 
        //var triangulatedCylinder = distMeshCylinder (origin, radiusCircle, heightCylinder);
    
        /*for (var m = 0; m<triangulatedCylinder.bars.length; m++) {
            const subArray = triangulatedCylinder.bars[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCylinder.finalPoints.length; m++) {
            const subArray = triangulatedCylinder.finalPoints[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCylinder.p.length; m++) {
            const subArray = triangulatedCylinder.p[m];
            console.log(subArray);
        }*/
    
        function distMeshCircle (origin, radius) {
            // pfix is a fixed set of points, which are given as input to the system.
            // The geometry is given as a distance function fd. This function returns the signed distance from each node location p to the closest boundary.
            // The (relative) desired edge length function h(x, y) is given as a function fh, which returns h for all input points.
            function fdCircle (points, radius) {
                var boundaryDistance = Math.sqrt(points[0]*points[0] + points[1]*points[1])-radius;
                return boundaryDistance;
            }
    
            function fhCircle (points) {
                var edgeLength = 2;
                return edgeLength;
            }
    
            var inputh0 = 2;
            
            var inputRadiusCircle = radius;
    
            var inputBboxCircle = [
                [-inputRadiusCircle,-inputRadiusCircle],
                [inputRadiusCircle,inputRadiusCircle]
            ];
            var inputPfixCircle = [
                [0,0], 
                [inputRadiusCircle,0], 
                [0,inputRadiusCircle],
                [-inputRadiusCircle,0], 
                [0,-inputRadiusCircle]
            ];
            var dptol=.001; 
            var ttol=.1; 
            var Fscale=1.2; 
            var deltat=.2; 
            var geps=.001*inputh0; 
            var deps=Math.sqrt(geps)*inputh0;
    
            // Creating an initial meshgrid
            const x_min = inputBboxCircle[0][0];
            const y_min = inputBboxCircle[0][1];
            const x_max = inputBboxCircle[1][0];
            const y_max = inputBboxCircle[1][1];
    
            const rows_x = Math.floor((x_max - x_min) / inputh0) + 1;
            const rows_y = Math.floor((y_max - y_min) / (inputh0 * Math.sqrt(3) / 2)) + 1;
            // rows_x and rows_y will contain the number of rows and columns needed in the provided bounding box dimensions based on the given value of h0. 
    
            const x = [];
            const y = [];
    
            for (let i = 0; i < rows_y; i++) {
                const row_x = [];
                const row_y = [];
            
                for (let j = 0; j < rows_x; j++) {
                    row_x.push(x_min + j * inputh0);
                    row_y.push(y_min + i * inputh0 * Math.sqrt(3) / 2);
                }
            
                x.push(row_x);
                y.push(row_y);
            }
            // From here, x and y will be an array of arrays. For the corresponding addresses, each elements in the x and y will be an array. 
            // The elements of the array within x and y will point to the corresponding coordinates. 
            // For instance, the array at the [0] index in both x and y will provide the x and y values corresponding to the coordinates at the first row of points. 
            
            // Shifting the x coordinates in the alternate rows by h0/2. 
            for (let i = 1; i < x.length; i += 2) {
                for (let j = 0; j < x[i].length; j++) {
                    x[i][j] += inputh0/2;
                }
            }
    
            // The .flat() functions converts the x and y from an array of arrays into a single array. 
            // All the elements are stacked in a series. 
            const flattenedX = x.flat();
            const flattenedY = y.flat();
            // console.log(flattenedX);
            //console.log(flattenedY.length);
    
            // Combine flattenedX and flattenedY column-wise
            // p developed from here will be an array of arrays, where each inside array will be a coordinate. 
            var p = flattenedX.map((value, index) => [value, flattenedY[index]]);
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // console.log(p.length);
            var filteredP = [];
    
            // Removing the points outside the boundary.
            // Only the interior points with negative distances (allowing a tolerance geps) are kept. 
            for (let i=0; i<p.length; i++) {
                const row = p[i];
    
                if (fdCircle(row, inputRadiusCircle) < geps) {
                    filteredP.push(row);
                }
            }
            // Again, the filteredP will also be same as P in terms of the type. It will contain the filtered points.  
            // console.log(filteredP.length);
    
            // Assigning filteredP back to p because all the further operations are done on p. 
            p = [...filteredP];
            // console.log(p.length);
            const r0 = [];
    
            // As per the paper, r0 is showing the probability to keep a point. So, it most probably will be just a value. 
            for (let i = 0; i < p.length; i++) {
                const row = p[i];
                const evaluation = fhCircle(row);
                const inverseSquare = 1 / evaluation ** 2;
                r0.push(inverseSquare);
            }
    
            // Generating a random array of the length of p.length. 
            const randValues = Array.from({ length: p.length }, () => Math.random());
            // console.log(randValues[0]);
    
            // Normalize the values in r0. 
            const normalizedR0 = r0.map(value => value / Math.max(...r0));
    
            // Selecting rows based on the condition
            const selectedRows = [];
    
            // I think, logically, this is doing nothing. I can probably remove this as well. What difference does the comparsion of a random value with its normalized value make? I am not sure with this. 
            for (let i = 0; i < p.length; i++) {
                if (randValues[i] < normalizedR0[i]) {
                selectedRows.push(p[i]);
                }
            }
    
            // Concatenating matrices. So, what I understand is, concatenatedP will also be an array of arrays, where each point coordinates are represented.
            const concatenatedP = inputPfixCircle.concat(selectedRows);
            p = [...concatenatedP];
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // Storing the length of p array in the variable N. 
            const N = p.length;
            //console.log(N);
    
            // Initializing pold with infinity
            let pold = Array.from({ length: p.length }, () => [Infinity, Infinity]);
            var numIter = 1;
            while (true) {
                // 3. Retriangulation by the Delaunay algorithm
                var displacement = 0;
                var difference = [];
                //Computing the elementwise difference between the p and pold. 
                for (var k=0; k<p.length; k++) {
                    difference[k] = [];
                    difference[k][0] = p[k][0] - pold[k][0];
                    difference[k][1] = p[k][1] - pold[k][1];
                }
    
                // Computing the square root of sum of squares of the difference and dividing by h0.
                var normalizedDifference = []; 
                for (var l=0; l<difference.length; l++) {
                    normalizedDifference[l] = Math.sqrt(difference[l][0]*difference[l][0] + difference[l][1]*difference[l][1])/inputh0;
                }
                // Finding the maximum value of the normalized difference from the array. 
                displacement = Math.max(...normalizedDifference);
                // console.log(displacement);
                //p.map((point, index) => Math.sqrt(point.reduce((sum, coord, axis) => sum + Math.pow(coord - pold[index][axis], 2), 0)) / h0);
            
                if (displacement > ttol) {
                    pold = p.map(point => [...point]); // Save current positions (deep copy). 
                    // We need to do this because if we directly assign pold = p, both will reference to the same object. So making changes in one would afftect the other as well. 
                    // console.log(pold[0]);
                    // So, since I am using the Delaunator library from the Javascript, the results are accordingly. 
                    // Delaunator.from(p).triangles will return an array. 
                    // Each consecutive three values in the triangles will show the indices of the points in p, which form a triangle in Delaunay.
                    // const trianglesIndices = Delaunator.from(p).triangles;
                    const trianglesIndices = new d3.Delaunay(p.flat()).triangles;
                    // console.log(trianglesIndices);
    
                    // Create an array of triangle points similar to what was needed in the original code, so that it can be processed accordingly.
                    const triangles = Array.from({ length: trianglesIndices.length/3 }, () => Array(3));
                    for (let i=0; i<trianglesIndices.length; i+=3){
                        for (let j=0; j<3; j++){
                            triangles[i/3][j] = trianglesIndices[i+j];
                        }
                    }
                    // console.log(triangles.length);
    
    
                    // Compute centroids of all the created triangles
                    const pmid = triangles.map(t => [
                        (p[t[0]][0] + p[t[1]][0] + p[t[2]][0]) / 3,
                        (p[t[0]][1] + p[t[1]][1] + p[t[2]][1]) / 3
                    ]);
                    // console.log(pmid[0]);
            
                    // Keep interior triangles. It filters those triangles whose centeroids are inside the boundary region. 
                    // Filter method will keep the triangles for which the condition is true. 
                    const interiorTriangles = triangles.filter((t, index) => fdCircle(pmid[index], inputRadiusCircle) < -geps); //Although the paper uses -geps here, I think it should be geps. 
                    // Looking into more detail, it does not make much difference geps or -geps). 
                    // console.log(interiorTriangles.length);
            
                    // Describe each bar by a unique pair of nodes. So basically from the selected nodes, you are creating the bars. 
                    var bars = interiorTriangles.flatMap(t => [
                        [t[0], t[1]],
                        [t[0], t[2]],
                        [t[1], t[2]]
                    ]);
                    // console.log(bars);
                    // To obtain a unique set of sorted bars from the bars saved in the previous step. 
                    const sortedBars = bars.map(bar => [...bar].sort((a, b) => a - b));
                    const uniqueBars = Array.from(new Set(sortedBars.map(JSON.stringify)), JSON.parse);
                    // Assigning uniqueBars back to the original bars, which is the variable used ahead in the code. 
                    bars = [...uniqueBars];
                    // console.log(bars);
                    // Graphical output of the current mesh (console.log used for demonstration)
                    //console.log("Iteration:", iterationCount);
                    //console.log("Triangles:", interiorTriangles);
                    //console.log("Bars:", uniqueBars);
            
                    //iterationCount++;
            
                    // Additional code for graphical output (adjust as needed)
                }
                // Now on, you will try to move the bar points based on the bar lengths and forces. 
                // So now, the indices stored in the bars need to be taken into the form of vector representing the difference between x and y components. 
                // List of bar vectors
                const barvec = bars.map(bar => [
                    p[bar[0]][0] - p[bar[1]][0],
                    p[bar[0]][1] - p[bar[1]][1]
                    ]);
                // console.log(barvec);
                // To my understanding, barvec will represent the x and y components of the vector representing the bars. 
    
    
                // Compute the length of each vector in the barvec.
                const L = barvec.map(vector => Math.sqrt(vector[0] ** 2 + vector[1] ** 2));
                // console.log(L);
    
                // Computing the midpoint of all the edges and then computing the fh function for these mid-point coordinates. 
                // So, hbars will be an array containing the function values of fh evaluated at the midpoints of the bars. 
                // Each element of hbars corresponds to a different bar.
                // Compute hbars
                const hbars = bars.map(bar => {
                    const midpoint = [
                    (p[bar[0]][0] + p[bar[1]][0]) / 2,
                    (p[bar[0]][1] + p[bar[1]][1]) / 2
                    ];
                    return fhCircle(midpoint);
                });
    
                // Computing the desired length of the bars
                const sumL2 = L.reduce((sum, value) => sum + value ** 2, 0);
                const sumHbars2 = hbars.reduce((sum, value) => sum + value ** 2, 0);
    
                const L0 = hbars.map((hbar) => hbar * Fscale * Math.sqrt(sumL2 / sumHbars2));
                // console.log(L0);
                // Computing the bar forces. It is being done by computing the difference between desired length and the computed length. 
                const F = L0.map((L0i, index) => Math.max(L0i - L[index], 0));
                // console.log (F);
                // Bar forces (x, y components)
                const Fvec = F.map((force, index) => [
                    force / L[index] * barvec[index][0],
                    force / L[index] * barvec[index][1]
                ]);
                // console.log(Fvec);
                //The resulting Fvec array contains arrays where each element represents the force vector (x and y components) for each bar.
                // console.log(Fvec);
                // Create a matrix representing the summation of forces at each node. 
    
                // var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                // The length of Ftot is same as that of the length of p (storing total number of nodes). This is to my understanding. 
    
                // Ftot has to be the sum of force vectors in Fvec, from all the bars meeting at a node.
                // A stretching force has a positive sign, and its direction is given by the two component vectors in bars. I think there is something wrong in this.  
                /*for (var z=0; z < p.length; z++) {
                    for (var a = 0; a < bars.length; a++) {
                        if (bars[a].includes(z)) {
                            Ftot[z][0] = Ftot[z][0] + Fvec[a][0];
                            Ftot[z][1] = Ftot[z][1] + Fvec[a][1];
                        }
                    }
                }*/
                // console.log(Ftot[0].length);
    
                var sparseInput1 = [];
                for (var i = 0; i < bars.length; i++) {
                    sparseInput1[i] = [];
                    sparseInput1[i][0] = bars[i][0];
                    sparseInput1[i][1] = bars[i][0];
                    sparseInput1[i][2] = bars[i][1];
                    sparseInput1[i][3] = bars[i][1];
                }
                var flatSparseInput1 = sparseInput1.flat();
    
                var originalArray = [0, 1, 0, 1];
                var lenF = F.length;
                var sparseInput2 = Array.from({ length: F.length }, () => [...originalArray]);
                var flatSparseInput2 = sparseInput2.flat();
    
                var sparseInput3 = [];
    
                for (var i = 0; i < Fvec.length; i++) {
                    sparseInput3[i] = [];
                    sparseInput3[i][0] = Fvec[i][0];
                    sparseInput3[i][1] = Fvec[i][1];
                    sparseInput3[i][2] = -1*Fvec[i][0];
                    sparseInput3[i][3] = -1*Fvec[i][1];
                }
                var flatSparseInput3 = sparseInput3.flat();
    
                var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                for (var i = 0; i < flatSparseInput1.length; i++) {
                    Ftot[flatSparseInput1[i]][flatSparseInput2[i]] += flatSparseInput3[i];
                }
    
    
                //bars.forEach((bar, index) => {
                //    Ftot[bar[0]][0] += Fvec[index][0];
                //    Ftot[bar[0]][1] += Fvec[index][1];
                //    Ftot[bar[1]][0] -= Fvec[index][0];
                //    Ftot[bar[1]][1] -= Fvec[index][1];
                //});
    
                // Setting up force components to be zero at the fixed nodes. 
                // Force = 0 at fixed points
                const numRows = inputPfixCircle.length; // Assuming pfix is a 2D array
                const numCols = Ftot[0].length; // Assuming Ftot is a 2D array
                for (let i = 0; i < numRows; i++) {
                    for (let j = 0; j < numCols; j++) {
                    Ftot[i][j] = 0;
                    }
                }
                //pfix.forEach(fixedPointIndex => {
                //    Ftot[fixedPointIndex][0] = 0;
                //    Ftot[fixedPointIndex][1] = 0;
                //});
    
                // Update node positions
                p = p.map((point, index) => [
                    point[0] + deltat * Ftot[index][0],
                    point[1] + deltat * Ftot[index][1]
                ]);
    
                // Finds points in p which are outside the d>0.
                // Resulting ix will be a boolean array, containing true if the distance > 0, false otherwise. 
                const d = p.map(point =>fdCircle(point, inputRadiusCircle));
                const ix = d.map(value => value > 0);
    
                // For the points where the value of d>0, you want to move them inwards. 
                // The first step is to find out the gradient of those points.
                // Then, based on the gradient, the points are moved. 
                const projectedPoints = ix.map((value, index) => {
                    if (value) {
                        // Here you make slight changes in the coordinate values at p. 
                        const perturbedX = [p[index][0] + deps, p[index][1]];
                        const perturbedY = [p[index][0], p[index][1] + deps];
                
                        // Here, you are trying to find out the gradient I guess. To see how much change should you do in the x and y coordinate. 
                        const dgradx = (fdCircle(perturbedX, inputRadiusCircle) - d[index]) / deps;
                        const dgrady = (fdCircle(perturbedY, inputRadiusCircle) - d[index]) / deps;
                        
                        // Here, you are making those changes in the x and y coordinate values. 
                        return [
                            p[index][0] - d[index] * dgradx,
                            p[index][1] - d[index] * dgrady
                        ];
                    } else {
                        return p[index];
                    }
                });
                // So the projected points here will be an array containing all the points in p, which also include the points which have been moved inwards. 
                // Update p with the projected points
                for (let i = 0; i < p.length; i++) {
                    p[i] = projectedPoints[i];
                }
    
                // Breaking criteria to come out of the loop. 
                // It checks if all the interior nodes move less than dptol. 
                // This is running infinitely long, I don't understand why. I am putting another break condition on the maximum number of iterations as 1000 based on https://github.com/ionhandshaker/distmesh/blob/master/README.md.
                const threshold = dptol;
                const logicalArray2 = d.map(element => element < -geps); // Create a logical array based on d for which the values are less than -geps. 
                var ftotLessThanGeps = [];
                for (let i = 0; i < logicalArray2.length; i++) {
                    if (logicalArray2[i]) {
                        ftotLessThanGeps.push(Ftot[i]);
                    }
                }
                var squaredNorms = [];
                for (let i = 0; i < ftotLessThanGeps.length; i++) {
                    squaredNorms[i] = Math.sqrt((ftotLessThanGeps[i][0]*deltat)**2 + (ftotLessThanGeps[i][1]*deltat)**2);
                }
                const maxSquaredNorm = Math.max(...squaredNorms);
                // const maxNorm = Math.sqrt(maxSquaredNorm);
                // console.log(maxSquaredNorm / h0);
                if (maxSquaredNorm / inputh0 < threshold) {
                    //return {p, bars}; 
                    break;
                } /*else if (numIter === 2000) {
                    return {p, bars}; 
                    break;
                } else {
                    numIter = numIter+1;
                }*/
                // Additional conditions and iterations may be added here
            }
            for (var i = 0;  i < p.length; i++) {
                p[i][2] = 0; 
                p[i][0] += origin[0];
                p[i][1] += origin[1];
                p[i][2] += origin[2];
            }
            var finalPoints = p;
            return {finalPoints, bars};
        }
    
        //var origin = [0,0,0];
        var radius = 10; 
        //var triangulatedCircle = distMeshCircle(origin, radius);
    
        /*for (var m = 0; m<triangulatedCircle.bars.length; m++) {
            const subArray = triangulatedCircle.bars[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCircle.finalPoints.length; m++) {
            const subArray = triangulatedCircle.finalPoints[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCylinder.cylinderPoints.length; m++) {
            const subArray = triangulatedCylinder.cylinderPoints[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedCircle.p.length; m++) {
            const subArray = triangulatedCircle.p[m];
            console.log(subArray);
        }*/
    
        function distMeshHollowCircle (origin, radius, innerRadius) {
            // pfix is a fixed set of points, which are given as input to the system.
            // The geometry is given as a distance function fd. This function returns the signed distance from each node location p to the closest boundary.
            // The (relative) desired edge length function h(x, y) is given as a function fh, which returns h for all input points.
            function fdHollowCircle (points, radius, innerRadius) {
                var boundaryDistance = Math.max((Math.sqrt(points[0]*points[0] + points[1]*points[1])-radius), -1*(Math.sqrt(points[0]*points[0] + points[1]*points[1])-innerRadius));
                return boundaryDistance;
            }
    
            function fhCircle (points) {
                var edgeLength = 2;
                return edgeLength;
            }
    
            var inputh0 = 2;
            
            var inputRadiusCircle = radius;
            var inputInternalRadiusCircle = innerRadius; 
    
            var inputBboxCircle = [
                [-inputRadiusCircle,-inputRadiusCircle],
                [inputRadiusCircle,inputRadiusCircle]
            ];
            var inputPfixCircle = [
                [0,0], 
                [inputRadiusCircle,0], 
                [0,inputRadiusCircle],
                [-inputRadiusCircle,0], 
                [0,-inputRadiusCircle]
            ];
            var dptol=.001; 
            var ttol=.1; 
            var Fscale=1.2; 
            var deltat=.2; 
            var geps=.001*inputh0; 
            var deps=Math.sqrt(geps)*inputh0;
    
            // Creating an initial meshgrid
            const x_min = inputBboxCircle[0][0];
            const y_min = inputBboxCircle[0][1];
            const x_max = inputBboxCircle[1][0];
            const y_max = inputBboxCircle[1][1];
    
            const rows_x = Math.floor((x_max - x_min) / inputh0) + 1;
            const rows_y = Math.floor((y_max - y_min) / (inputh0 * Math.sqrt(3) / 2)) + 1;
            // rows_x and rows_y will contain the number of rows and columns needed in the provided bounding box dimensions based on the given value of h0. 
    
            const x = [];
            const y = [];
    
            for (let i = 0; i < rows_y; i++) {
                const row_x = [];
                const row_y = [];
            
                for (let j = 0; j < rows_x; j++) {
                    row_x.push(x_min + j * inputh0);
                    row_y.push(y_min + i * inputh0 * Math.sqrt(3) / 2);
                }
            
                x.push(row_x);
                y.push(row_y);
            }
            // From here, x and y will be an array of arrays. For the corresponding addresses, each elements in the x and y will be an array. 
            // The elements of the array within x and y will point to the corresponding coordinates. 
            // For instance, the array at the [0] index in both x and y will provide the x and y values corresponding to the coordinates at the first row of points. 
            
            // Shifting the x coordinates in the alternate rows by h0/2. 
            for (let i = 1; i < x.length; i += 2) {
                for (let j = 0; j < x[i].length; j++) {
                    x[i][j] += inputh0/2;
                }
            }
    
            // The .flat() functions converts the x and y from an array of arrays into a single array. 
            // All the elements are stacked in a series. 
            const flattenedX = x.flat();
            const flattenedY = y.flat();
            // console.log(flattenedX);
            //console.log(flattenedY.length);
    
            // Combine flattenedX and flattenedY column-wise
            // p developed from here will be an array of arrays, where each inside array will be a coordinate. 
            var p = flattenedX.map((value, index) => [value, flattenedY[index]]);
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // console.log(p.length);
            var filteredP = [];
    
            // Removing the points outside the boundary.
            // Only the interior points with negative distances (allowing a tolerance geps) are kept. 
            for (let i=0; i<p.length; i++) {
                const row = p[i];
    
                if (fdHollowCircle(row, inputRadiusCircle, inputInternalRadiusCircle) < geps) {
                    filteredP.push(row);
                }
            }
            // Again, the filteredP will also be same as P in terms of the type. It will contain the filtered points.  
            // console.log(filteredP.length);
    
            // Assigning filteredP back to p because all the further operations are done on p. 
            p = [...filteredP];
            // console.log(p.length);
            const r0 = [];
    
            // As per the paper, r0 is showing the probability to keep a point. So, it most probably will be just a value. 
            for (let i = 0; i < p.length; i++) {
                const row = p[i];
                const evaluation = fhCircle(row);
                const inverseSquare = 1 / evaluation ** 2;
                r0.push(inverseSquare);
            }
    
            // Generating a random array of the length of p.length. 
            const randValues = Array.from({ length: p.length }, () => Math.random());
            // console.log(randValues[0]);
    
            // Normalize the values in r0. 
            const normalizedR0 = r0.map(value => value / Math.max(...r0));
    
            // Selecting rows based on the condition
            const selectedRows = [];
    
            // I think, logically, this is doing nothing. I can probably remove this as well. What difference does the comparsion of a random value with its normalized value make? I am not sure with this. 
            for (let i = 0; i < p.length; i++) {
                if (randValues[i] < normalizedR0[i]) {
                selectedRows.push(p[i]);
                }
            }
    
            // Concatenating matrices. So, what I understand is, concatenatedP will also be an array of arrays, where each point coordinates are represented.
            const concatenatedP = inputPfixCircle.concat(selectedRows);
            p = [...concatenatedP];
            /*for (var o = 0; o< p.length; o++) {
                console.log(p[o]); 
            }*/
            
            // Storing the length of p array in the variable N. 
            const N = p.length;
            //console.log(N);
    
            // Initializing pold with infinity
            let pold = Array.from({ length: p.length }, () => [Infinity, Infinity]);
            var numIter = 1;
            while (true) {
                // 3. Retriangulation by the Delaunay algorithm
                var displacement = 0;
                var difference = [];
                //Computing the elementwise difference between the p and pold. 
                for (var k=0; k<p.length; k++) {
                    difference[k] = [];
                    difference[k][0] = p[k][0] - pold[k][0];
                    difference[k][1] = p[k][1] - pold[k][1];
                }
    
                // Computing the square root of sum of squares of the difference and dividing by h0.
                var normalizedDifference = []; 
                for (var l=0; l<difference.length; l++) {
                    normalizedDifference[l] = Math.sqrt(difference[l][0]*difference[l][0] + difference[l][1]*difference[l][1])/inputh0;
                }
                // Finding the maximum value of the normalized difference from the array. 
                displacement = Math.max(...normalizedDifference);
                // console.log(displacement);
                //p.map((point, index) => Math.sqrt(point.reduce((sum, coord, axis) => sum + Math.pow(coord - pold[index][axis], 2), 0)) / h0);
            
                if (displacement > ttol) {
                    pold = p.map(point => [...point]); // Save current positions (deep copy). 
                    // We need to do this because if we directly assign pold = p, both will reference to the same object. So making changes in one would afftect the other as well. 
                    // console.log(pold[0]);
                    // So, since I am using the Delaunator library from the Javascript, the results are accordingly. 
                    // Delaunator.from(p).triangles will return an array. 
                    // Each consecutive three values in the triangles will show the indices of the points in p, which form a triangle in Delaunay.
                    // const trianglesIndices = Delaunator.from(p).triangles;
                    const trianglesIndices = new d3.Delaunay(p.flat()).triangles;
                    // console.log(trianglesIndices);
    
                    // Create an array of triangle points similar to what was needed in the original code, so that it can be processed accordingly.
                    const triangles = Array.from({ length: trianglesIndices.length/3 }, () => Array(3));
                    for (let i=0; i<trianglesIndices.length; i+=3){
                        for (let j=0; j<3; j++){
                            triangles[i/3][j] = trianglesIndices[i+j];
                        }
                    }
                    // console.log(triangles.length);
    
    
                    // Compute centroids of all the created triangles
                    const pmid = triangles.map(t => [
                        (p[t[0]][0] + p[t[1]][0] + p[t[2]][0]) / 3,
                        (p[t[0]][1] + p[t[1]][1] + p[t[2]][1]) / 3
                    ]);
                    // console.log(pmid[0]);
            
                    // Keep interior triangles. It filters those triangles whose centeroids are inside the boundary region. 
                    // Filter method will keep the triangles for which the condition is true. 
                    var interiorTriangles = triangles.filter((t, index) => fdHollowCircle(pmid[index], inputRadiusCircle, inputInternalRadiusCircle) < -geps); //Although the paper uses -geps here, I think it should be geps. 
                    // Looking into more detail, it does not make much difference geps or -geps). 
                    // console.log(interiorTriangles.length);
            
                    // Describe each bar by a unique pair of nodes. So basically from the selected nodes, you are creating the bars. 
                    var bars = interiorTriangles.flatMap(t => [
                        [t[0], t[1]],
                        [t[0], t[2]],
                        [t[1], t[2]]
                    ]);
                    // console.log(bars);
                    // To obtain a unique set of sorted bars from the bars saved in the previous step. 
                    const sortedBars = bars.map(bar => [...bar].sort((a, b) => a - b));
                    const uniqueBars = Array.from(new Set(sortedBars.map(JSON.stringify)), JSON.parse);
                    // Assigning uniqueBars back to the original bars, which is the variable used ahead in the code. 
                    bars = [...uniqueBars];
                    // console.log(bars);
                    // Graphical output of the current mesh (console.log used for demonstration)
                    //console.log("Iteration:", iterationCount);
                    //console.log("Triangles:", interiorTriangles);
                    //console.log("Bars:", uniqueBars);
            
                    //iterationCount++;
            
                    // Additional code for graphical output (adjust as needed)
                }
                // Now on, you will try to move the bar points based on the bar lengths and forces. 
                // So now, the indices stored in the bars need to be taken into the form of vector representing the difference between x and y components. 
                // List of bar vectors
                const barvec = bars.map(bar => [
                    p[bar[0]][0] - p[bar[1]][0],
                    p[bar[0]][1] - p[bar[1]][1]
                    ]);
                // console.log(barvec);
                // To my understanding, barvec will represent the x and y components of the vector representing the bars. 
    
    
                // Compute the length of each vector in the barvec.
                const L = barvec.map(vector => Math.sqrt(vector[0] ** 2 + vector[1] ** 2));
                // console.log(L);
    
                // Computing the midpoint of all the edges and then computing the fh function for these mid-point coordinates. 
                // So, hbars will be an array containing the function values of fh evaluated at the midpoints of the bars. 
                // Each element of hbars corresponds to a different bar.
                // Compute hbars
                const hbars = bars.map(bar => {
                    const midpoint = [
                    (p[bar[0]][0] + p[bar[1]][0]) / 2,
                    (p[bar[0]][1] + p[bar[1]][1]) / 2
                    ];
                    return fhCircle(midpoint);
                });
    
                // Computing the desired length of the bars
                const sumL2 = L.reduce((sum, value) => sum + value ** 2, 0);
                const sumHbars2 = hbars.reduce((sum, value) => sum + value ** 2, 0);
    
                const L0 = hbars.map((hbar) => hbar * Fscale * Math.sqrt(sumL2 / sumHbars2));
                // console.log(L0);
                // Computing the bar forces. It is being done by computing the difference between desired length and the computed length. 
                const F = L0.map((L0i, index) => Math.max(L0i - L[index], 0));
                // console.log (F);
                // Bar forces (x, y components)
                const Fvec = F.map((force, index) => [
                    force / L[index] * barvec[index][0],
                    force / L[index] * barvec[index][1]
                ]);
                // console.log(Fvec);
                //The resulting Fvec array contains arrays where each element represents the force vector (x and y components) for each bar.
                // console.log(Fvec);
                // Create a matrix representing the summation of forces at each node. 
    
                // var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                // The length of Ftot is same as that of the length of p (storing total number of nodes). This is to my understanding. 
    
                // Ftot has to be the sum of force vectors in Fvec, from all the bars meeting at a node.
                // A stretching force has a positive sign, and its direction is given by the two component vectors in bars. I think there is something wrong in this.  
                /*for (var z=0; z < p.length; z++) {
                    for (var a = 0; a < bars.length; a++) {
                        if (bars[a].includes(z)) {
                            Ftot[z][0] = Ftot[z][0] + Fvec[a][0];
                            Ftot[z][1] = Ftot[z][1] + Fvec[a][1];
                        }
                    }
                }*/
                // console.log(Ftot[0].length);
    
                var sparseInput1 = [];
                for (var i = 0; i < bars.length; i++) {
                    sparseInput1[i] = [];
                    sparseInput1[i][0] = bars[i][0];
                    sparseInput1[i][1] = bars[i][0];
                    sparseInput1[i][2] = bars[i][1];
                    sparseInput1[i][3] = bars[i][1];
                }
                var flatSparseInput1 = sparseInput1.flat();
    
                var originalArray = [0, 1, 0, 1];
                var lenF = F.length;
                var sparseInput2 = Array.from({ length: F.length }, () => [...originalArray]);
                var flatSparseInput2 = sparseInput2.flat();
    
                var sparseInput3 = [];
    
                for (var i = 0; i < Fvec.length; i++) {
                    sparseInput3[i] = [];
                    sparseInput3[i][0] = Fvec[i][0];
                    sparseInput3[i][1] = Fvec[i][1];
                    sparseInput3[i][2] = -1*Fvec[i][0];
                    sparseInput3[i][3] = -1*Fvec[i][1];
                }
                var flatSparseInput3 = sparseInput3.flat();
    
                var Ftot = new Array(N).fill(null).map(() => [0, 0]); // creates an array of size N, with all elements as [0,0].
                for (var i = 0; i < flatSparseInput1.length; i++) {
                    Ftot[flatSparseInput1[i]][flatSparseInput2[i]] += flatSparseInput3[i];
                }
    
    
                //bars.forEach((bar, index) => {
                //    Ftot[bar[0]][0] += Fvec[index][0];
                //    Ftot[bar[0]][1] += Fvec[index][1];
                //    Ftot[bar[1]][0] -= Fvec[index][0];
                //    Ftot[bar[1]][1] -= Fvec[index][1];
                //});
    
                // Setting up force components to be zero at the fixed nodes. 
                // Force = 0 at fixed points
                const numRows = inputPfixCircle.length; // Assuming pfix is a 2D array
                const numCols = Ftot[0].length; // Assuming Ftot is a 2D array
                for (let i = 0; i < numRows; i++) {
                    for (let j = 0; j < numCols; j++) {
                    Ftot[i][j] = 0;
                    }
                }
                //pfix.forEach(fixedPointIndex => {
                //    Ftot[fixedPointIndex][0] = 0;
                //    Ftot[fixedPointIndex][1] = 0;
                //});
    
                // Update node positions
                p = p.map((point, index) => [
                    point[0] + deltat * Ftot[index][0],
                    point[1] + deltat * Ftot[index][1]
                ]);
    
                // Finds points in p which are outside the d>0.
                // Resulting ix will be a boolean array, containing true if the distance > 0, false otherwise. 
                const d = p.map(point => fdHollowCircle(point, inputRadiusCircle, inputInternalRadiusCircle));
                const ix = d.map(value => value > 0);
    
                // For the points where the value of d>0, you want to move them inwards. 
                // The first step is to find out the gradient of those points.
                // Then, based on the gradient, the points are moved. 
                const projectedPoints = ix.map((value, index) => {
                    if (value) {
                        // Here you make slight changes in the coordinate values at p. 
                        const perturbedX = [p[index][0] + deps, p[index][1]];
                        const perturbedY = [p[index][0], p[index][1] + deps];
                
                        // Here, you are trying to find out the gradient I guess. To see how much change should you do in the x and y coordinate. 
                        const dgradx = (fdHollowCircle(perturbedX, inputRadiusCircle, inputInternalRadiusCircle) - d[index]) / deps;
                        const dgrady = (fdHollowCircle(perturbedY, inputRadiusCircle, inputInternalRadiusCircle) - d[index]) / deps;
                        
                        // Here, you are making those changes in the x and y coordinate values. 
                        return [
                            p[index][0] - d[index] * dgradx,
                            p[index][1] - d[index] * dgrady
                        ];
                    } else {
                        return p[index];
                    }
                });
                // So the projected points here will be an array containing all the points in p, which also include the points which have been moved inwards. 
                // Update p with the projected points
                for (let i = 0; i < p.length; i++) {
                    p[i] = projectedPoints[i];
                }
    
                // Breaking criteria to come out of the loop. 
                // It checks if all the interior nodes move less than dptol. 
                // This is running infinitely long, I don't understand why. I am putting another break condition on the maximum number of iterations as 1000 based on https://github.com/ionhandshaker/distmesh/blob/master/README.md.
                const threshold = dptol;
                const logicalArray2 = d.map(element => element < -geps); // Create a logical array based on d for which the values are less than -geps. 
                var ftotLessThanGeps = [];
                for (let i = 0; i < logicalArray2.length; i++) {
                    if (logicalArray2[i]) {
                        ftotLessThanGeps.push(Ftot[i]);
                    }
                }
                var squaredNorms = [];
                for (let i = 0; i < ftotLessThanGeps.length; i++) {
                    squaredNorms[i] = Math.sqrt((ftotLessThanGeps[i][0]*deltat)**2 + (ftotLessThanGeps[i][1]*deltat)**2);
                }
                const maxSquaredNorm = Math.max(...squaredNorms);
                // const maxNorm = Math.sqrt(maxSquaredNorm);
                // console.log(maxSquaredNorm / h0);
                if (maxSquaredNorm / inputh0 < threshold) {
                    //return {p, bars}; 
                    break;
                } /*else if (numIter === 2000) {
                    return {p, bars}; 
                    break;
                } else {
                    numIter = numIter+1;
                }*/
                // Additional conditions and iterations may be added here
            }
            for (var i = 0;  i < p.length; i++) {
                p[i][2] = 0; 
                p[i][0] += origin[0];
                p[i][1] += origin[1];
                p[i][2] += origin[2];
            }
            var finalPoints = p;
            return {finalPoints, bars, interiorTriangles};
        }
    
        // var origin = [0,0,0];
        var radius = 10;
        var innerRadius = 3;
        // var triangulatedHollowCircle = distMeshHollowCircle(origin, radius, innerRadius);
    
        /*for (var m = 0; m<triangulatedHollowCircle.bars.length; m++) {
            const subArray = triangulatedHollowCircle.bars[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedHollowCircle.finalPoints.length; m++) {
            const subArray = triangulatedHollowCircle.finalPoints[m];
            console.log(subArray);
        }*/
        /*for (var m = 0; m<triangulatedHollowCircle.p.length; m++) {
            const subArray = triangulatedHollowCircle.p[m];
            console.log(subArray);
        }*/
    
        function getSurfaceEquation(polyhedralSurface) {
            // Returns the a,b,c,d for the plane equation of the form ax+by+cz+d = 0. 
            // Also returns the coefficients for the unit normal vector n1, n2, n3. 
            let coefficients = [];
            let lineCoefficient1 = [];
            let lineCoefficient2 = [];
            let normalVector = [];
            lineCoefficient1[0] = polyhedralSurface[1][0]-polyhedralSurface[0][0];
            lineCoefficient1[1] = polyhedralSurface[1][1]-polyhedralSurface[0][1];
            lineCoefficient1[2] = polyhedralSurface[1][2]-polyhedralSurface[0][2];
            lineCoefficient2[0] = polyhedralSurface[3][0]-polyhedralSurface[0][0];
            lineCoefficient2[1] = polyhedralSurface[3][1]-polyhedralSurface[0][1];
            lineCoefficient2[2] = polyhedralSurface[3][2]-polyhedralSurface[0][2];
            normalVector[0] = lineCoefficient1[1]*lineCoefficient2[2] - lineCoefficient2[1]*lineCoefficient1[2];
            normalVector[1] = -1 * (lineCoefficient1[0]*lineCoefficient2[2] - lineCoefficient2[0]*lineCoefficient1[2]);
            normalVector[2] = lineCoefficient1[0]*lineCoefficient2[1] - lineCoefficient2[0]*lineCoefficient1[1];
            let d = 0 - (normalVector[0]*polyhedralSurface[0][0] + normalVector[1]*polyhedralSurface[0][1] + normalVector[2]*polyhedralSurface[0][2]);
            let n1 = normalVector[0]/Math.sqrt(normalVector[0]**2 + normalVector[1]**2 + normalVector[2]**2);
            let n2 = normalVector[1]/Math.sqrt(normalVector[0]**2 + normalVector[1]**2 + normalVector[2]**2);
            let n3 = normalVector[2]/Math.sqrt(normalVector[0]**2 + normalVector[1]**2 + normalVector[2]**2);
            let Jd = 0 - (n1*polyhedralSurface[0][0] + n2*polyhedralSurface[0][1] + n3*polyhedralSurface[0][2]);
            return {normalVector0: normalVector[0], normalVector1: normalVector[1], normalVector2: normalVector[2], d, n1, n2, n3, Jd};
        }
    
        function findSignedDistance (surface, edgeCoordinates) {
            // This function computes the Signed distance between the two ends of the edge and the surface. 
            let signedDistance = [];
            // Using the formula used by Thomas Moller in his paper.
            signedDistance[0] = surface.d + (surface.normalVector0*edgeCoordinates[0][0] + surface.normalVector1*edgeCoordinates[0][1] + surface.normalVector2*edgeCoordinates[0][2]);
            signedDistance[1] = surface.d + (surface.normalVector0*edgeCoordinates[1][0] + surface.normalVector1*edgeCoordinates[1][1] + surface.normalVector2*edgeCoordinates[1][2]);
            return signedDistance;
        }
    
        function checkIfIntersectionOccurs(surface, edgeCoordinates) {
            // This function finds if the intersection occurs between the surface and the bar. 
            let signedDistance = findSignedDistance (surface, edgeCoordinates);
            if (signedDistance[0]==0 && signedDistance[1]==0) {
                return "coplanar";
            }
            else if ((signedDistance[0]==0 && signedDistance[1]!=0) || (signedDistance[0]!=0 && signedDistance[1]==0)) {
                return "one point on the surface";
            } 
            else if ((signedDistance[0]<0 && signedDistance[1]>0) || (signedDistance[0]>0 && signedDistance[1]<0)) {
                return true;
            }
            else {
                return false;
            }
        }
    
        function getLineEquation (edgeCoordinates) {
            // This function returns the equation of the line based on the coordinates of the two end points. 
            // It considers the line to be in the form of x = x0 + ta, y = y0 +tb, z = z0 + tc. 
            // Here, a, b, and c are the coefficients of the vector in the direction parallel to the line. 
            // The function returns the value of x0, y0, z0, a, b, and c. 
            let a = edgeCoordinates[1][0] - edgeCoordinates[0][0];
            let b = edgeCoordinates[1][1] - edgeCoordinates[0][1];
            let c = edgeCoordinates[1][2] - edgeCoordinates[0][2];
            let x0 = edgeCoordinates[0][0];
            let y0 = edgeCoordinates[0][1];
            let z0 = edgeCoordinates[0][2];
            return {a, b, c, x0, y0, z0};
        }
    
        function findIntersectionCoordinates(surfaceEquation, edgeCoordinates, intersectionCheck, polyhedralSurface){
    
            // The function to find the coordinates of intersection between the input surface and the edge having coordinates of the two ends. 
            // The value of the variable intersectionCheck is important to find the intersection coordinates accordingly. 
            // The value of the parameter t is found by solving the line equation in parametric form and the plane equation. 
            // The value of t is being found using t = - (D + Ax0 + By0 + Cz0)/(Aa + Bb + Cc)
            let coordinates = [];
            if (intersectionCheck==true) {
                // If the two points are on the opposite side, find the intersection point using t. 
                let lineEquation = getLineEquation(edgeCoordinates);
                let t = -1*((surfaceEquation.d + surfaceEquation.normalVector0*lineEquation.x0 + surfaceEquation.normalVector1*lineEquation.y0 + surfaceEquation.normalVector2*lineEquation.z0)/(surfaceEquation.normalVector0*lineEquation.a + surfaceEquation.normalVector1*lineEquation.b + surfaceEquation.normalVector2*lineEquation.c));
                coordinates[0] = []
                coordinates[0][0] = lineEquation.x0 + t * lineEquation.a ;
                coordinates[0][1] = lineEquation.y0 + t * lineEquation.b ;
                coordinates[0][2] = lineEquation.z0 + t * lineEquation.c ; 
            }
            else if (intersectionCheck == "one point on the surface") {
                // If one point is on the surface, then return the coordinates of that point itself.
                coordinates[0] = []
                let signedDistance = findSignedDistance (surfaceEquation, edgeCoordinates);
                if (signedDistance[0]==0 && signedDistance[1]!=0) {
                    coordinates[0][0] = edgeCoordinates[0][0];
                    coordinates[0][1] = edgeCoordinates[0][1];
                    coordinates[0][2] = edgeCoordinates[0][2];
                }
                else {
                    coordinates[0][0] = edgeCoordinates[1][0];
                    coordinates[0][1] = edgeCoordinates[1][1];
                    coordinates[0][2] = edgeCoordinates[1][2];
                }
            }
            else if (intersectionCheck == "coplanar") {
                // I am not sure about this case. If the line lies on the surface, it means every point is intersecting. So, taking just the end points might not be correct. 
                // Will have a look at it later. 
                // If both points are there, this in itself means the whole line is there.
                // I am not sure at this point if this will really make a difference or not. 
                coordinates[0] = [];
                coordinates[1] = [];
                coordinates[0][0] = edgeCoordinates[0][0];
                coordinates[0][1] = edgeCoordinates[0][1];
                coordinates[0][2] = edgeCoordinates[0][2];
                coordinates[1][0] = edgeCoordinates[1][0];
                coordinates[1][1] = edgeCoordinates[1][1];
                coordinates[1][2] = edgeCoordinates[1][2];
                // An important point for this case will be finding the coordinates of intersection of the bar and the edges of the polyhedral surface. 
                for (let i=0; i<polyhedralSurface.length-1; i++) {
                    let line1 = getLineEquation([polyhedralSurface[i], polyhedralSurface[i+1]]);
                    const x01 = line1.x0;
                    const y01 = line1.y0;
                    const z01 = line1.z0;
                    const A1 = line1.a;
                    const B1 = line1.b;
                    const C1 = line1.c;
                    let line2 = getLineEquation(edgeCoordinates);
                    let A2 = line2.a;
                    let B2 = line2.b;
                    let C2 = line2.c;
                    let x02 = line2.x0;
                    let y02 = line2.y0;
                    let z02 = line2.z0;
                    if (((B2*A1-A2*B1)!=0) || ((B1*C2-B2*C1)!=0) || ((A1*C2-A2*C1)!=0)){
                        // If any of the above is zero, the line will be parallel to that line. But still the intersection may exist. 
                            if ((B2*A1-A2*B1)!=0) {
                                var t1 = (B2*(x02-x01) - A2*(y02-y01))/(B2*A1-A2*B1);
                                var t2 = (B1*(x02-x01) - A1*(y02-y01))/(B2*A1-A2*B1);
                            }
                            else if ((B1*C2-B2*C1)!=0) {
                                var t1 = (C2*(y02-y01) - B2*(z02-z01))/(B1*C2-B2*C1);
                                var t2 = (C1*(y02-y01) - B1*(z02-z01))/(B1*C2-B2*C1);
                            }
                            else if ((A1*C2-A2*C1)!=0) {
                                var t1 = (C2*(x02-x01) - A2*(z02-z01))/(A1*C2-A2*C1);
                                var t2 = (C1*(x02-x01) - A1*(z02-z01))/(A1*C2-A2*C1);
                            }
                            // Computing the coordinates of the intersection point
                            // The solution is found only by solving the x and y equations. So, to check if the solution actually exists, we need to check if the z equation also satisfies the identified values. 
            
                            if (((z01+C1*t1).toFixed(4)===(z02+C2*t2).toFixed(4)) && ((x01+A1*t1).toFixed(4)===(x02+A2*t2).toFixed(4)) && ((y01+B1*t1).toFixed(4)===(y02+B2*t2).toFixed(4)) && t1 >= 0 && t1<=1 && t2 >= 0 && t2 <= 1) {
                                let x_intersect = x01 + A1 * t1;
                                let y_intersect = y01 + B1 * t1;
                                let z_intersect = z01 + C1 * t1;
                                coordinates.push([x_intersect, y_intersect, z_intersect]);
                            }
                        }    
                }
            }
            return coordinates;
        }
    
        function getTrianglePlaneEquation(triangulatedSurface, trianglePointsIndex) {
            // This function will return the equation of plane for the triangle, same as the earlier one for the cube surface. 
            let coefficients = [];
            let lineCoefficient1 = [];
            let lineCoefficient2 = [];
            let normalVector = [];
            lineCoefficient1[0] = triangulatedSurface.finalPoints[trianglePointsIndex[1]][0]-triangulatedSurface.finalPoints[trianglePointsIndex[0]][0];
            lineCoefficient1[1] = triangulatedSurface.finalPoints[trianglePointsIndex[1]][1]-triangulatedSurface.finalPoints[trianglePointsIndex[0]][1];
            lineCoefficient1[2] = triangulatedSurface.finalPoints[trianglePointsIndex[1]][2]-triangulatedSurface.finalPoints[trianglePointsIndex[0]][2];
            lineCoefficient2[0] = triangulatedSurface.finalPoints[trianglePointsIndex[2]][0]-triangulatedSurface.finalPoints[trianglePointsIndex[0]][0];
            lineCoefficient2[1] = triangulatedSurface.finalPoints[trianglePointsIndex[2]][1]-triangulatedSurface.finalPoints[trianglePointsIndex[0]][1];
            lineCoefficient2[2] = triangulatedSurface.finalPoints[trianglePointsIndex[2]][2]-triangulatedSurface.finalPoints[trianglePointsIndex[0]][2];
            normalVector[0] = lineCoefficient1[1]*lineCoefficient2[2] - lineCoefficient2[1]*lineCoefficient1[2];
            normalVector[1] = -1 * (lineCoefficient1[0]*lineCoefficient2[2] - lineCoefficient2[0]*lineCoefficient1[2]);
            normalVector[2] = lineCoefficient1[0]*lineCoefficient2[1] - lineCoefficient2[0]*lineCoefficient1[1];
            let d = 0 - (normalVector[0]*triangulatedSurface.finalPoints[trianglePointsIndex[0]][0] + normalVector[1]*triangulatedSurface.finalPoints[trianglePointsIndex[0]][1] + normalVector[2]*triangulatedSurface.finalPoints[trianglePointsIndex[0]][2]);
            let n1 = normalVector[0]/Math.sqrt(normalVector[0]**2 + normalVector[1]**2 + normalVector[2]**2);
            let n2 = normalVector[1]/Math.sqrt(normalVector[0]**2 + normalVector[1]**2 + normalVector[2]**2);
            let n3 = normalVector[2]/Math.sqrt(normalVector[0]**2 + normalVector[1]**2 + normalVector[2]**2);
            let Jd = 0 - (n1*triangulatedSurface.finalPoints[trianglePointsIndex[0]][0] + n2*triangulatedSurface.finalPoints[trianglePointsIndex[0]][0] + n3*triangulatedSurface.finalPoints[trianglePointsIndex[0]][0]);
            return {normalVector0: normalVector[0], normalVector1: normalVector[1], normalVector2: normalVector[2], d, n1, n2, n3, Jd};
    
        }
    
        function checkIfLiesWithinLimit(surface, coordinates, shape) {
            //console.log(surface);
            //console.log(coordinates);
            if(shape == "TriangleEdgetoCubeSurface") {
                // This function is to check if a point lies within the boundaries of a Surface. 
                // The surface can be anything, a rectangle or a triangle. 
                // Ray casting algorithm is used to this end.
                const x01 = coordinates[0];
                const y01 = coordinates[1];
                const z01 = coordinates[2];
                let inside = false;
                // Getting the equation of a line in any direction in the plane. 
                // I am considering the ray to be passed along the line formed by joining the first two coordinates of the surface.
                // Equation of the line is of the form Ai + Bj + Ck
                const A1 = surface[1][0] - surface[0][0];
                const B1 = surface[1][1] - surface[0][1];
                const C1 = surface[1][2] - surface[0][2];
                // I am assuming the value of t as 100000 here, as I think it should be fine.
                //const x = x0 + A*t;
                //const y = y0 + B*t;
                //const z = z0 + C*t;
                // The code is somehow not working properly for the points that lie on the nodes or the edges. 
                // So I am giving explicit conditions for them. Otherwise, the logic looks fine to me.
                // First, to check if the coordinate lies on the nodes.  
                for (let i=0; i<surface.length; i++) {
                    if ((coordinates[0]==surface[i][0]) && (coordinates[1]==surface[i][1]) && (coordinates[2]==surface[i][2])){
                        return true;
                    }
                }
                // Next, to check if the coordinate lies on the edges. 
                // A point can be considered on a line if a solution for t exists for that point in the parametric form of the line. 
                // If the point lies on the line, check if it is within the bounds of the line. 
                // If not within the bounds, return false, so that there is no further check. 
                for (let i=0; i<surface.length-1; i++) {
                    let Aline1 = x01-surface[i][0];
                    let Bline1 = y01-surface[i][1];
                    let Cline1 = z01-surface[i][2];
                    let Aline2 = surface[i+1][0]-x01;
                    let Bline2 = surface[i+1][1]-y01;
                    let Cline2 = surface[i+1][2]-z01;
                    let line1Xline2 = [];
                    line1Xline2[0] = Bline1*Cline2-Bline2*Cline1;
                    line1Xline2[1] = -1* (Aline1*Cline2 - Aline2*Cline1);
                    line1Xline2[2] = Aline1*Bline2 - Aline2*Bline1;
                    if((line1Xline2[0]==0 || line1Xline2[0]==-0) && 
                        (line1Xline2[1]==0 || line1Xline2[1]==-0) && 
                        (line1Xline2[2]==0 || line1Xline2[2]==-0) &&
                        (Math.min(surface[i][0],surface[i+1][0])<=x01 && x01<=Math.max(surface[i][0],surface[i+1][0])) && 
                        (Math.min(surface[i][1],surface[i+1][1])<=y01 && y01<=Math.max(surface[i][1],surface[i+1][1])) && 
                        (Math.min(surface[i][2],surface[i+1][2])<=z01 && z01<=Math.max(surface[i][2],surface[i+1][2]))) {
                            // console.log(1);
                            return true;
                        }
                    }
    
    
                var intersectionCount = 0; 
                for (let i=0; i<surface.length-1; i++) {
                    let line = getLineEquation([surface[i], surface[i+1]]);
                    let A2 = line.a;
                    let B2 = line.b;
                    let C2 = line.c;
                    let x02 = line.x0;
                    let y02 = line.y0;
                    let z02 = line.z0;
                    if (((B2*A1-A2*B1)!=0) || ((B1*C2-B2*C1)!=0) || ((A1*C2-A2*C1)!=0)){
                    // If any of the above is zero, the line will be parallel to that line. But still the intersection may exist. 
                        if ((B2*A1-A2*B1)!=0) {
                            var t1 = (B2*(x02-x01) - A2*(y02-y01))/(B2*A1-A2*B1);
                            var t2 = (B1*(x02-x01) - A1*(y02-y01))/(B2*A1-A2*B1);
                        }
                        else if ((B1*C2-B2*C1)!=0) {
                            var t1 = (C2*(y02-y01) - B2*(z02-z01))/(B1*C2-B2*C1);
                            var t2 = (C1*(y02-y01) - B1*(z02-z01))/(B1*C2-B2*C1);
                        }
                        else if ((A1*C2-A2*C1)!=0) {
                            var t1 = (C2*(x02-x01) - A2*(z02-z01))/(A1*C2-A2*C1);
                            var t2 = (C1*(x02-x01) - A1*(z02-z01))/(A1*C2-A2*C1);
                        }
                        // Computing the coordinates of the intersection point
                        // The solution is found only by solving the x and y equations. So, to check if the solution actually exists, we need to check if the z equation also satisfies the identified values. 
    
                        if ((((z01+C1*t1).toFixed(4)===(z02+C2*t2).toFixed(4)) || ((z01+C1*t1).toFixed(4)==0 && (z02+C2*t2).toFixed(4)==0)) 
                            && (((x01+A1*t1).toFixed(4)===(x02+A2*t2).toFixed(4)) || ((x01+A1*t1).toFixed(4)==0 && (x02+A2*t2).toFixed(4)==0))
                            && ((y01+B1*t1).toFixed(4)===(y02+B2*t2).toFixed(4) || ((y01+B1*t1).toFixed(4)==0 && (y02+B2*t2).toFixed(4)==0)) 
                            && t1 >= 0 && t2 >= 0 && t2 <= 1) {
                            let x_intersect = x02 + A2 * t2;
                            let y_intersect = y02 + B2 * t2;
                            let z_intersect = z02 + C2 * t2;
                            // Now we need to check if the intersection point lies within the bounds of the line. 
                            // If it lies within the line, the intersectionCount is increased by 1. 
                            if ((Math.min(surface[i][0], surface[i+1][0])<=x_intersect && x_intersect<=Math.max(surface[i][0], surface[i+1][0])) && (Math.min(surface[i][1], surface[i+1][1])<=y_intersect && y_intersect<=Math.max(surface[i][1], surface[i+1][1])) && (Math.min(surface[i][2], surface[i+1][2])<=z_intersect && z_intersect<=Math.max(surface[i][2], surface[i+1][2]))) {
                                intersectionCount = intersectionCount + 1;
                            }
                        }
                    }
                }   
            }
    
            else if(shape == "CubeEdgetoTriangleSurface") {
                // This function is to check if a point lies within the boundaries of a Surface. 
                // The surface can be anything, a rectangle or a triangle. 
                // Ray casting algorithm is used to this end.
                const x01 = coordinates[0];
                const y01 = coordinates[1];
                const z01 = coordinates[2];
                let inside = false;
                // Getting the equation of a line in any direction in the plane. 
                // I am considering the ray to be passed along the line formed by joining the first two coordinates of the surface.
                // Equation of the line is of the form Ai + Bj + Ck
                const A1 = surface[1][0] - surface[0][0];
                const B1 = surface[1][1] - surface[0][1];
                const C1 = surface[1][2] - surface[0][2];
                // I am assuming the value of t as 100000 here, as I think it should be fine.
                //const x = x0 + A*t;
                //const y = y0 + B*t;
                //const z = z0 + C*t;
                // The code is somehow not working properly for the points that lie on the nodes or the edges. 
                // So I am giving explicit conditions for them. Otherwise, the logic looks fine to me.
                // First, to check if the coordinate lies on the nodes.  
                for (let i=0; i<surface.length; i++) {
                    // Check if the coordinate lies on the nodes of the triangle.
                    if ((coordinates[0]==surface[i][0]) && (coordinates[1]==surface[i][1]) && (coordinates[2]==surface[i][2])){
                        return true;
                    }
                }
                // Next, to check if the coordinate lies on the edges. 
                // For this, first check if the lines formed by connecting the intersection coordinate and the two ends of the edge have a zero cross product. 
                // Then check if the coordinate lies within the limits of the two ends of the edge.  
                // If not within the bounds, return false, so that there is no further check. 
                for (let i=0; i<surface.length-1; i++) {
                    let Aline1 = x01-surface[i][0];
                    let Bline1 = y01-surface[i][1];
                    let Cline1 = z01-surface[i][2];
                    let Aline2 = surface[i+1][0]-x01;
                    let Bline2 = surface[i+1][1]-y01;
                    let Cline2 = surface[i+1][2]-z01;
                    let line1Xline2 = [];
                    line1Xline2[0] = Bline1*Cline2-Bline2*Cline1;
                    line1Xline2[1] = -1* (Aline1*Cline2 - Aline2*Cline1);
                    line1Xline2[2] = Aline1*Bline2 - Aline2*Bline1;
                    if((line1Xline2[0]==0 || line1Xline2[0]==-0) && 
                        (line1Xline2[1]==0 || line1Xline2[1]==-0) && 
                        (line1Xline2[2]==0 || line1Xline2[2]==-0) &&
                        (Math.min(surface[i][0],surface[i+1][0])<=x01 && x01<=Math.max(surface[i][0],surface[i+1][0])) && 
                        (Math.min(surface[i][1],surface[i+1][1])<=y01 && y01<=Math.max(surface[i][1],surface[i+1][1])) && 
                        (Math.min(surface[i][2],surface[i+1][2])<=z01 && z01<=Math.max(surface[i][2],surface[i+1][2]))) {
                            // console.log(1);
                            return true;
                        }
                    }
    
                // Next, if the intersection coordinate lie somewhere on the surface. 
                var intersectionCount = 0; 
                for (let i=0; i<surface.length-1; i++) {
                    let line = getLineEquation([surface[i], surface[i+1]]);
                    let A2 = line.a;
                    let B2 = line.b;
                    let C2 = line.c;
                    let x02 = line.x0;
                    let y02 = line.y0;
                    let z02 = line.z0;
                    if (((B2*A1-A2*B1)!=0) || ((B1*C2-B2*C1)!=0) || ((A1*C2-A2*C1)!=0)){
                    // If any of the above is zero, the line will be parallel to that line. But still the intersection may exist. 
                        if ((B2*A1-A2*B1)!=0) {
                            var t1 = (B2*(x02-x01) - A2*(y02-y01))/(B2*A1-A2*B1);
                            var t2 = (B1*(x02-x01) - A1*(y02-y01))/(B2*A1-A2*B1);
                        }
                        else if ((B1*C2-B2*C1)!=0) {
                            var t1 = (C2*(y02-y01) - B2*(z02-z01))/(B1*C2-B2*C1);
                            var t2 = (C1*(y02-y01) - B1*(z02-z01))/(B1*C2-B2*C1);
                        }
                        else if ((A1*C2-A2*C1)!=0) {
                            var t1 = (C2*(x02-x01) - A2*(z02-z01))/(A1*C2-A2*C1);
                            var t2 = (C1*(x02-x01) - A1*(z02-z01))/(A1*C2-A2*C1);
                        }
                        // Computing the coordinates of the intersection point
                        // The solution is found only by solving the x and y equations. So, to check if the solution actually exists, we need to check if the z equation also satisfies the identified values. 
    
                        if ((((z01+C1*t1).toFixed(4)===(z02+C2*t2).toFixed(4)) || ((z01+C1*t1).toFixed(4)==0 && (z02+C2*t2).toFixed(4)==0)) 
                            && (((x01+A1*t1).toFixed(4)===(x02+A2*t2).toFixed(4)) || ((x01+A1*t1).toFixed(4)==0 && (x02+A2*t2).toFixed(4)==0))
                            && ((y01+B1*t1).toFixed(4)===(y02+B2*t2).toFixed(4) || ((y01+B1*t1).toFixed(4)==0 && (y02+B2*t2).toFixed(4)==0)) 
                            && t1 >= 0 && t2 >= 0 && t2 <= 1) {
                            let x_intersect = x02 + A2 * t2;
                            let y_intersect = y02 + B2 * t2;
                            let z_intersect = z02 + C2 * t2;
                            // Now we need to check if the intersection point lies within the bounds of the line. 
                            // If it lies within the line, the intersectionCount is increased by 1. 
                            if ((Math.min(surface[i][0], surface[i+1][0])<=x_intersect && x_intersect<=Math.max(surface[i][0], surface[i+1][0])) && (Math.min(surface[i][1], surface[i+1][1])<=y_intersect && y_intersect<=Math.max(surface[i][1], surface[i+1][1])) && (Math.min(surface[i][2], surface[i+1][2])<=z_intersect && z_intersect<=Math.max(surface[i][2], surface[i+1][2]))) {
                                intersectionCount = intersectionCount + 1;
                            }
                        }
                    }
                }   
            }
            // Need to add the case of CSpace checking in this function. 
            else if (shape == "polygon") {
                // You need to assume that if the shape is polygon, then its function is being called to check the point and C-Space
                // Another thing is, here you just need to work in 2D. Because whatever we are checking is all on the plans. 
                // Defining the line parallel to X-axis from the coordinate under consideration
                // We need to add a check if the point lies on any of the edges. In that case, we directly increase the intersection count. 
                let intersectionCoordinates = [];
                let line1 = {x0 : coordinates[0],
                            y0 : coordinates[1],
                            a : 1,
                            b : 0};
                let A1 = line1.a;
                let B1 = line1.b;
                let x01 = line1.x0;
                let y01 = line1.y0;
                // iterate over all the edges to check for the intersection. 
                for (let i=0; i<surface.edges.length; i++) {
                    let edge = surface.edges[i];
                    let line2 = getLineEquation([surface.coordinates[edge[0]], surface.coordinates[edge[1]]]);
                    if (checkIfLiesInLine(coordinates, ([surface.coordinates[edge[0]], surface.coordinates[edge[1]]])).status == true) {
                        return true;          
                    }
                    else {
                        let A2 = line2.a;
                        let B2 = line2.b;
                        let x02 = line2.x0;
                        let y02 = line2.y0;
                        if ((B2*A1-A2*B1)!=0) {
                            var t1 = (B2*(x02-x01) - A2*(y02-y01))/(B2*A1-A2*B1);
                            var t2 = (B1*(x02-x01) - A1*(y02-y01))/(B2*A1-A2*B1);
                            if (t1 >= 0 && t2 >= 0 && t2 <= 1) {
                                // So, when you check that B2*A1-A2*B1!=0, this means that the lines will always intersect. 
                                // Checking t1>=0 means that we are considering only one direction. 
                                // For t2<=1, we check that the coordinates of intersection will lie within the bounds of the line segment. 
                                let x_intersect = x02 + A2 * t2;
                                let y_intersect = y02 + B2 * t2;
                                intersectionCoordinates.push([x_intersect, y_intersect]);
                            }
                        }
                    }
                }
                // Now you need to remove duplicates from the list
                // Use a Set to store unique string representations of inner arrays
                let uniqueIntersectionCoordinates = new Set();
                // Iterate through each inner array
                for (let i = 0; i < intersectionCoordinates.length; i++) {
                    // Convert the inner array to a string for comparison
                    let arrayAsString = JSON.stringify(intersectionCoordinates[i]);
                    // Check if the string representation is already in the Set
                    if (!uniqueIntersectionCoordinates.has(arrayAsString)) {
                        // If not, add it to the Set
                        uniqueIntersectionCoordinates.add(arrayAsString);
                    }
                }
                // Convert the Set back to an array of arrays
                let uniqueIntersectionCoordinatesList = Array.from(uniqueIntersectionCoordinates, JSON.parse);
                intersectionCount = uniqueIntersectionCoordinatesList.length;
            }
            // If the intersection is done in even numbers, the point is outside, else inside.
            if ((intersectionCount%2==0) || (intersectionCount==0)) {
                return false;
            }
            else {
                return true;
            }
        }
    
        function getPlaneDistanceFromPoint (surface, point) {
            let plane1 = getSurfaceEquation(surface);
            let distance = Math.abs(plane1.normalVector0*point[0]+plane1.normalVector1*point[1]+plane1.normalVector2*point[2]+plane1.d)/Math.sqrt(plane1.normalVector0**2+plane1.normalVector1**2+plane1.normalVector2**2);
            return distance;
        }
    
        function getLineDistanceFromPoint (line, point){
            let APx1 = point[0]-line[0][0];
            let APy1 = point[1]-line[0][1];
            let APz1 = point[2]-line[0][2];
            let dx1 = line[1][0] - line[0][0];
            let dy1 = line[1][1] - line[0][1];
            let dz1 = line[1][2] - line[0][2];
            let CP1 = [];
            CP1[0] = (APy1*dz1-dy1*APz1);
            CP1[1] = -1*(APx1*dz1-dx1*APz1);
            CP1[2] = (APx1*dy1-dx1*APy1);
            let D1 = Math.abs(Math.sqrt(CP1[0]**2+CP1[1]**2+CP1[2]**2)/Math.sqrt(dx1**2+dy1**2+dz1**2));
            return D1;
        }
    
        function detectIntersection (triangulatedSurface, polyhedralSurface) {
            // triangulatedSurface should contain the points, bars, and the set of points forming each triangle. 
            // polyhedralSurface should be the typical polyhedralSurface representation containing the coordinates of corners of each face.
            // We will do a two-way detection. 
            // First, we will detect the intersection between the triangle edges and polyhedral surface. 
            // Next, we will detect the intersection between the edges of the polyhedral surface and the triangle. 
            // The coordinates of intersection will be found for both separately.
            // It will return a list containing the coordinates of the point of intersection. 
            var intersectionCoordinates = [];
            var intersectionCoordinatesEdges = [];
            // Creating a separate function for detecting the intersection between the triangle edges and cube surfaces first. 
            function detectTriangleEdgetoCubeSurface (triangulatedSurface, polyhedralSurface) {
                // Defining the equation of the plane for each surface of the polyhedral surface.
                // planeEquation list will store the coefficients a,b,c,d for the equation of the plane.
                let shape = "TriangleEdgetoCubeSurface" ;
                var planeEquation = [];
                for (let s = 0; s < polyhedralSurface.length; s++) {
                    planeEquation.push(getSurfaceEquation(polyhedralSurface[s]));
                }
                for (var b=0; b<triangulatedSurface.bars.length; b++) {
                    // We can just look at the bars and detect the intersection between the bars and polyhedral surface.  
                    // Every time this loop runs, the value of barCoordinates gets updated. At a time, it will store just two coordinates. 
                    var barCoordinates = [];
                    for (let j=0; j<2; j++) {
                        barCoordinates.push(triangulatedSurface.finalPoints[triangulatedSurface.bars[b][j]]);
                    }
                    // We need to find the surface equation for each face of the polyhedral surface. 
                    // We also need to find the equation of the line. 
                    // Now, we check the bar for intersection with each surface. 
                    for (let pl=0; pl<planeEquation.length; pl++) {
                        // First, check if the intersection occurs between the bar and the surface. 
                        let intersectionCheck = checkIfIntersectionOccurs(planeEquation[pl], barCoordinates);
                        if ((intersectionCheck == true) || (intersectionCheck == "coplanar") || (intersectionCheck == "one point on the surface")) {
                            // If intersection occurs, find the intersection coordinates. 
                            let coordinatesList = findIntersectionCoordinates(planeEquation[pl], barCoordinates, intersectionCheck, polyhedralSurface[pl]);
                            for (let i=0; i<coordinatesList.length; i++) {
                                // Need to check if the identified intersection point lies within the limit of the cube surface.
                                let coordinate = checkIfLiesWithinLimit(polyhedralSurface[pl], coordinatesList[i], shape);
                                if (coordinate == true) {
                                    // Although, I do not think this makes sense, but I am also adding a check to make sure that the coordinate lies within the bounds of the intersecting line. 
                                    let A = barCoordinates[1][0]-barCoordinates[0][0];
                                    let B = barCoordinates[1][1]-barCoordinates[0][1];
                                    let C = barCoordinates[1][2]-barCoordinates[0][2];
                                    if (A!=0){
                                        let t = (coordinatesList[i][0]-barCoordinates[0][0])/A;
                                        if (t>=0 && t<=1) {
                                            intersectionCoordinates.push(coordinatesList[i]);
                                            intersectionCoordinatesEdges.push(barCoordinates);
                                        }
                                    }
                                    else if (B!=0){
                                        let t = (coordinatesList[i][1]-barCoordinates[0][1])/B;
                                        if (t>=0 && t<=1) {
                                            intersectionCoordinates.push(coordinatesList[i]);
                                            intersectionCoordinatesEdges.push(barCoordinates);
                                        }
                                    }
                                    else if (C!=0){
                                        let t = (coordinatesList[i][2]-barCoordinates[0][2])/C;
                                        if (t>=0 && t<=1) {
                                            intersectionCoordinates.push(coordinatesList[i]);
                                            intersectionCoordinatesEdges.push(barCoordinates);
                                        }
                                    }
                                }
                            }  
                        }
                    }   
                }
            }
            // Creating a separate function to detect the intersection between the triangle surface and cube edges.
            function detectCubeEdgetoTriangleSurface (triangulatedSurface, polyhedralSurface){
                // For every set of coordinates in the polyhedral surface, performing the intersection checking. 
                // Check is performed for all the triangles for each set of coordinates in the polyhedral surface.
                //You need to make sure that the function is capable of handling the data when it comes in WKT triangulated format.
                
                // Getting each polyhedral surface from the set of surfaces. 
                let shape = "CubeEdgetoTriangleSurface";
                for (let s=0; s<polyhedralSurface.length; s++){
    
                    // Iterate over the interior triangles to select each triangle respectively. 
                    for (let j=0; j<triangulatedSurface.interiorTriangles.length; j++){
    
                        // Getting the triangle plane equation for the selected triangle. 
                        var trianglePlaneEquation = getTrianglePlaneEquation(triangulatedSurface, triangulatedSurface.interiorTriangles[j]);
                        let triangleSurface = [triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][0]],triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][1]],triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][2]],triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][0]]];
    
                        // Getting the coordinates of each line in the polyhedral surface. 
                        for (let c=0; c<polyhedralSurface[s].length-1; c++) {
    
                            // check if the intersection occurs between the triangle plane equation and each line in the polyhedral surface.
                            let intersectionCheck = checkIfIntersectionOccurs(trianglePlaneEquation, [polyhedralSurface[s][c],polyhedralSurface[s][c+1]]);
    
                            if (intersectionCheck == true || intersectionCheck == "coplanar" || intersectionCheck == "one point on the surface") {
                                let coordinatesList = findIntersectionCoordinates(trianglePlaneEquation, [polyhedralSurface[s][c],polyhedralSurface[s][c+1]], intersectionCheck, triangleSurface);
                                // let triangleCoordinates = [triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][0]],triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][1]],triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][2]], triangulatedSurface.finalPoints[triangulatedSurface.interiorTriangles[j][0]]];
    
                                for (let i=0; i<coordinatesList.length; i++) {
                                    let coordinate = checkIfLiesWithinLimit(triangleSurface,coordinatesList[i],shape);
                                    if (coordinate == true) {
                                        // Although, I do not think this makes sense, but I am also adding a check to make sure that the coordinate lies within the bounds of the intersecting line.
                                        let A = polyhedralSurface[s][c+1][0]-polyhedralSurface[s][c][0];
                                        let B = polyhedralSurface[s][c+1][1]-polyhedralSurface[s][c][1];
                                        let C = polyhedralSurface[s][c+1][2]-polyhedralSurface[s][c][2];
                                        if (A!=0){
                                            let t = (coordinatesList[i][0]-polyhedralSurface[s][c][0])/A;
                                            if (t>=0 && t<=1) {
                                                intersectionCoordinates.push(coordinatesList[i]);
                                                intersectionCoordinatesEdges.push([polyhedralSurface[s][c],polyhedralSurface[s][c+1]]);
                                            }
                                        }
                                        else if (B!=0){
                                            let t = (coordinatesList[i][1]-polyhedralSurface[s][c][1])/B;
                                            if (t>=0 && t<=1) {
                                                intersectionCoordinates.push(coordinatesList[i]);
                                                intersectionCoordinatesEdges.push([polyhedralSurface[s][c],polyhedralSurface[s][c+1]]);
                                            }
                                        }
                                        else if (C!=0){
                                            let t = (coordinatesList[i][2]-polyhedralSurface[s][c][2])/C;
                                            if (t>=0 && t<=1) {
                                                intersectionCoordinates.push(coordinatesList[i]);
                                                intersectionCoordinatesEdges.push([polyhedralSurface[s][c],polyhedralSurface[s][c+1]]);
                                            }
                                        }
                                    }   
                                }
                            }   
                        }
                    }
                }
            }
            detectTriangleEdgetoCubeSurface (triangulatedSurface, polyhedralSurface);
            detectCubeEdgetoTriangleSurface (triangulatedSurface, polyhedralSurface);
            return {intersectionCoordinates, intersectionCoordinatesEdges};
            // Detection of the intersection points being inside the boundary is also remaining.
        }
    
        function checkIfPointLiesOnPlane (point, plane) {
            if ((parseFloat((plane.normalVector0*point[0]+plane.normalVector1*point[1]+plane.normalVector2*point[2]+plane.d).toFixed(4))===0) ||
            (parseFloat((plane.normalVector0*point[0]+plane.normalVector1*point[1]+plane.normalVector2*point[2]+plane.d).toFixed(4))===-0)) {
                return true;
            }
        }
    
        function checkIfLinesAreparallel (edge1, edge2) {
            let A1 = edge1[1][0] - edge1[0][0];
            let B1 = edge1[1][1] - edge1[0][1];
            let C1 = edge1[1][2] - edge1[0][2];
            let A2 = edge2[1][0] - edge2[0][0];
            let B2 = edge2[1][1] - edge2[0][1];
            let C2 = edge2[1][2] - edge2[0][2];
            let crossProduct = [
                B1 * C2 - B2 * C1,
                A2 * C1 - A1 * C2,
                A1 * B2 - A2 * B1
            ];
            if ((crossProduct[0]==0 || crossProduct[0]==-0) && (crossProduct[1]==0 || crossProduct[1]==-0) && (crossProduct[2]==0 || crossProduct[2]==-0)) {
                return true;
            }
            else {
                return false;
            }  
        }
    
        function checkIfSurfacesAreparallel (surface1, surface2) {
            let surfaceEqn1 = getSurfaceEquation(surface1);
            let surfaceEqn2 = getSurfaceEquation(surface2);
            let A1 = surfaceEqn1.normalVector0;
            let B1 = surfaceEqn1.normalVector1;
            let C1 = surfaceEqn1.normalVector2;
            let A2 = surfaceEqn2.normalVector0;
            let B2 = surfaceEqn2.normalVector1;
            let C2 = surfaceEqn2.normalVector2;
            let crossProduct = [
                B1 * C2 - B2 * C1,
                A2 * C1 - A1 * C2,
                A1 * B2 - A2 * B1
            ];
            if ((crossProduct[0]==0 || crossProduct[0]==-0) && (crossProduct[1]==0 || crossProduct[1]==-0) && (crossProduct[2]==0 || crossProduct[2]==-0)) {
                return true;
            }
            else {
                return false;
            }
        }
    
        function getCloserCoordinates (edges, origin) {
            let closestEnds = [];
            let closestEndsIndex = [];
            let dist1 = Math.sqrt((edges[0][0][0]-origin[0])**2 + (edges[0][0][1]-origin[1])**2 + (edges[0][0][2]-origin[2])**2);
            let dist2 = Math.sqrt((edges[0][1][0]-origin[0])**2 + (edges[0][1][1]-origin[1])**2 + (edges[0][1][2]-origin[2])**2);
            let dist3 = Math.sqrt((edges[1][0][0]-origin[0])**2 + (edges[1][0][1]-origin[1])**2 + (edges[1][0][2]-origin[2])**2);
            let dist4 = Math.sqrt((edges[1][1][0]-origin[0])**2 + (edges[1][1][1]-origin[1])**2 + (edges[1][1][2]-origin[2])**2);
    
            let distances = [dist1, dist2, dist3, dist4];
            let lowest1 = Infinity;
            let lowest2 = Infinity;
    
            // Getting the closest two ends. 
            for (const num of distances) {
                if (num < lowest1) {
                lowest2 = lowest1;
                lowest1 = num;
                } else if (num < lowest2) {
                lowest2 = num;
                }
            }
            
            // Getting the indices of the closest distances. 
            for (let i=0; i<distances.length; i++) {
                if (distances[i]==lowest1 || distances[i]==lowest2) {
                    closestEndsIndex.push(i);
                }
            }
    
            if (closestEndsIndex.includes(0)){
                closestEnds.push(edges[0][0]);
            }
            if (closestEndsIndex.includes(1)){
                closestEnds.push(edges[0][0]);
            }
            if (closestEndsIndex.includes(2)){
                closestEnds.push(edges[1][0]);
            }
            if (closestEndsIndex.includes(3)){
                closestEnds.push(edges[1][1]);
            }
            return closestEnds;
        }
    
        /*var polyhedralSurface = [[[-2, 2, 4],[-2, 6, 4],[2, 6, 4], [2, 2, 4], [-2, 2, 4]],
                                [[-2, 2, 4], [-2, 2, 0], [2, 2, 0], [2, 2, 4], [-2, 2, 4]],
                                [[-2, 6, 4], [-2, 6, 0], [2, 6, 0], [2, 6, 4], [-2, 6, 4]],
                                [[-2, 2, 4], [-2, 2, 0], [-2, 6, 0], [-2, 6, 4], [-2, 2, 4]],
                                [[2, 2, 4], [2, 6, 4], [2, 6, 0], [2, 2, 0], [2, 2, 4]]];*///[[-3.51251, 3.558703, 4.000167 ],[-5, 3.558703, 4.000167],[-5, 6, 4.000167],[-3.51251, 6, 4.000167],[-3.51251, 3.558703, 4.000167]]];//,[[2,0,0],[4,0,0],[4,0,3],[2,0,3],[2,0,0]],[[2,0,3],[2,3,3],[4,3,3],[4,0,3],[2,0,3]], [[2,0,0],[2,0,3],[2,3,3],[2,3,0],[2,0,0]]]//,[[2,3,5],[2,3,0],[2,0,0],[2,0,5],[2,3,5]],[[4,3,5],[4,3,0],[4,0,0],[4,0,5],[4,3,5]]]//[[2,0,4],[2,3,4],[4,3,4],[4,0,4],[2,0,4]],[[2,0,0],[4,0,0],[4,3,0],[2,3,0],[2,0,0]],[[2,3,4],[4,3,4],[4,3,0],[2,3,0],[2,3,4]]];[[1,-6,3],[2,-6,3],[2,6,3],[1,6,3],[1,-6,3]],
        // var triangulatedSurface = distMeshCone (origin, radiusCircle, heightCone);
        /*for (let i=0; i<triangulatedSurface.finalPoints.length; i++) {
            console.log(triangulatedSurface.finalPoints[i]);
        }*/
        // console.log(triangulatedSurface.interiorTriangles);
    
    
        /*for (let i=0; i<triangulatedSurface.bars.length; i++) {
            console.log(triangulatedSurface.bars[i]);
        }*/
        //var coordinatesOfIntersection = detectIntersection (triangulatedSurface, polyhedralSurface);
        //console.log(coordinatesOfIntersection.intersectionCoordinates);
        // onsole.log(coordinatesOfIntersection.intersectionCoordinatesEdges);
    
        function checkIfLiesInLine (coordinate, edgeCoordinates) {
            if (edgeCoordinates[0].length == 3) {
                let a1 = coordinate[0]-edgeCoordinates[0][0];
                let b1 = coordinate[1]-edgeCoordinates[0][1];
                let c1 = coordinate[2]-edgeCoordinates[0][2];
                let a2 = edgeCoordinates[1][0] - coordinate[0];
                let b2 = edgeCoordinates[1][1] - coordinate[1];
                let c2 = edgeCoordinates[1][2] - coordinate[2];
                let corssProduct = [];
                corssProduct[0] = b1*c2-b2*c1;
                corssProduct[1] = -1*(a1*c2-a2*c1);
                corssProduct[2] = a1*b2 - a2*b1;
                if (((corssProduct[0]==0) || (corssProduct[0]==-0)) &&  ((corssProduct[1]==0) || (corssProduct[1]==-0)) && ((corssProduct[2]==0) || (corssProduct[2]==-0))) {
                    if (((Math.min(edgeCoordinates[0][0],edgeCoordinates[1][0])) <= coordinate[0]) && (coordinate[0]<=(Math.max(edgeCoordinates[0][0],edgeCoordinates[1][0]))) &&
                    ((Math.min(edgeCoordinates[0][1],edgeCoordinates[1][1])) <= coordinate[1]) && (coordinate[1]<=(Math.max(edgeCoordinates[0][1],edgeCoordinates[1][1]))) &&
                    ((Math.min(edgeCoordinates[0][2],edgeCoordinates[1][2])) <= coordinate[2]) && (coordinate[2]<=(Math.max(edgeCoordinates[0][2],edgeCoordinates[1][2])))){
                        let status = true;
                        let edge = edgeCoordinates;
                        return {status, edge};
                    }
                }
                else {
                    let status = false;
                    let edge = edgeCoordinates;
                    return {status, edge};
                }
            }
            else if (edgeCoordinates[0].length == 2) {
                // edgeCoordinates[0].push(0);
                // edgeCoordinates[1].push(0);
                let a1 = coordinate[0]-edgeCoordinates[0][0];
                let b1 = coordinate[1]-edgeCoordinates[0][1];
                let c1 = 0;
                let a2 = edgeCoordinates[1][0] - coordinate[0];
                let b2 = edgeCoordinates[1][1] - coordinate[1];
                let c2 = 0;
                let corssProduct = [];
                corssProduct[0] = b1*c2-b2*c1;
                corssProduct[1] = -1*(a1*c2-a2*c1);
                corssProduct[2] = a1*b2 - a2*b1;
                if (((corssProduct[0]==0) || (corssProduct[0]==-0)) &&  ((corssProduct[1]==0) || (corssProduct[1]==-0)) && ((corssProduct[2]==0) || (corssProduct[2]==-0))) {
                    if (((Math.min(edgeCoordinates[0][0],edgeCoordinates[1][0])) <= coordinate[0]) && (coordinate[0]<=(Math.max(edgeCoordinates[0][0],edgeCoordinates[1][0]))) &&
                    ((Math.min(edgeCoordinates[0][1],edgeCoordinates[1][1])) <= coordinate[1]) && (coordinate[1]<=(Math.max(edgeCoordinates[0][1],edgeCoordinates[1][1])))){
                        let status = true;
                        let edge = edgeCoordinates;
                        return {status, edge};
                    }
                    else {
                        let status = false;
                        let edge = edgeCoordinates;
                        return {status, edge};
                    }
                }
                else {
                    let status = false;
                    let edge = edgeCoordinates;
                    return {status, edge};
                }
            }    
        }
    
        function getExtremeCoordinates (coordinateOfIntersection, coordinateOfIntersectionEdges, triangulatedSurface, polyhedralSurface, heightOfRotatingObject, origin, radiusCircle) {
            let zList = [];
            // Getting the unique coordinates of intersection. 
            let intersectionCoordinates = Array.from(new Set(coordinateOfIntersection.map(JSON.stringify)), JSON.parse);
            // The intersectionCoordinatesAtTop will store the coordinates of intersection. 
            let intersectionCoordinatesAtTop = [];
            // The intersectionCoordinateEdges will store the corresponding edge coordinates. 
            let intersectionCoordinateEdges = [];
            let intersectionCoordinatesAtTopUnique = [];
            let intersectionCoordinateEdgesUnique = [];
            let extremeCoordinates = [];
            let extremeCoordinateEdges = [];
            let extremeCoordinateSurfaces = [];
            // Getting the maximum heigth of all the coordinates of the polyhedralsurface. 
            for (let i=0; i<polyhedralSurface.length; i++) {
                for (let j=0; j<polyhedralSurface[i].length; j++) {
                    zList.push(polyhedralSurface[i][j][2]);
                } 
            }
            let heightPolyhedralSurface = Math.max(...zList);
            // If the height of crane boom is greater than height of polyhedral surface, the extreme points will lie on the height of the polyhedral surface. 
            if (heightOfRotatingObject >= heightPolyhedralSurface) {
                for (let i=0; i<polyhedralSurface.length; i++) {
                    for (let j=0; j<polyhedralSurface[i].length-1; j++) {
                        // Check if the edge is the top edge. We check if the height of both the ends of the edge is same and is equal to the maximum height. 
                        if ((polyhedralSurface[i][j][2]===polyhedralSurface[i][j+1][2]) && (polyhedralSurface[i][j][2]===heightPolyhedralSurface)) {
                            let line = getLineEquation(polyhedralSurface[i][j], polyhedralSurface[i][j+1]);
                            // Find which intersection coordinates lie on the top edges. If they lie, add them to the intersectionCoordinatesAtTop
                            for (let k=0; k<intersectionCoordinates.length; k++) {
                                let pointOnLine = checkIfLiesInLine (intersectionCoordinates[k], [polyhedralSurface[i][j],polyhedralSurface[i][j+1]]);
                                if (pointOnLine.status == true) {
                                    // If the point lies on an edge, we store the point and the edge coordinates in the defined arrays. 
                                    // To avoid any problem due to numerical errors, I am rounding off them to 4 digits and doing parsefloat to convert them back to numbers. 
                                    intersectionCoordinatesAtTop.push(intersectionCoordinates[k].map((element)=>parseFloat(element.toFixed(4))));
                                    intersectionCoordinateEdges.push (pointOnLine.edge);
                                }
                            }
                        }
                    }
                }
                // Need to remove duplicates from the arrays prepared above. This is because the same edge might occur multiple times as it might be common in different cube faces. 
                // Another thing is edge coordinates may alse be reverse but same. This also has to be taken care off.
                for (let i=0; i<intersectionCoordinatesAtTop.length; i++) {
                    if (i==0){
                        intersectionCoordinatesAtTopUnique.push(intersectionCoordinatesAtTop[i]);
                        intersectionCoordinateEdgesUnique.push(intersectionCoordinateEdges[i]);
                    }
                    if (i!=0){
                        let count = 0;
                        for (let j=0; j<i; j++) {
                            if ((JSON.stringify(intersectionCoordinatesAtTop[i])==JSON.stringify(intersectionCoordinatesAtTop[j])) && 
                            ((JSON.stringify(intersectionCoordinateEdges[i])==JSON.stringify(intersectionCoordinateEdges[j])) ||
                            (JSON.stringify(intersectionCoordinateEdges[i])==JSON.stringify([intersectionCoordinateEdges[j][1],intersectionCoordinateEdges[j][0]])))) {
                                count+=1;
                            }
                        }
                        if (count==0) {
                            intersectionCoordinatesAtTopUnique.push(intersectionCoordinatesAtTop[i]);
                            intersectionCoordinateEdgesUnique.push(intersectionCoordinateEdges[i]);
                        }
                    } 
                }
                // Generally, the intersectionCoordinatesAtTopUnique will have two points. 
                // However, in some cases it can have one point (if the cube just touces the triangular surface). 
                // In some cases, it can also have four points (if the two edges of the bounding box intersect the cube, each at two different places). 
                // In case of four intersection points, we probably need to find the two which will be the governing ones.
                // An important thing in this case will be if there are four coordinates but only two lines, this means they lie on two parallel lines of the bounding box.  
                // So, first we need to find out the two parallel lines. 
                // The coordinates on the line which is orthogonally closer to the origin will govern the coordinates. 
                if (intersectionCoordinatesAtTopUnique.length==4) {
                    let uniqueEdges = [];
                    // let uniqueCoordinates = [];
                    // uniqueCoordinates.push(intersectionCoordinatesAtTopUnique[0]);
                    // Getting all the unique edges from all the edges.
                    // Unique edges is an array which will include only the unique edges which have intersection points. 
                    uniqueEdges.push(intersectionCoordinateEdgesUnique[0]);
                    for (let i=1; i<intersectionCoordinatesAtTopUnique.length; i++) {
                        let count = 0; 
                        for (let j=0; j<i; j++) {
                            if (JSON.stringify(intersectionCoordinateEdgesUnique[i])===JSON.stringify(intersectionCoordinateEdgesUnique[j])){
                                count+=1;
                            }
                        }
                        if (count==0){
                            uniqueEdges.push(intersectionCoordinateEdgesUnique[i]);
                        }
                    }
                    // If the length of the unique edges is 2, this means these two are the parallel edges. 
                    // Finding the closer edge. Because that will be the governing edge. 
                    // Using the formula |APXd|/|d| to find out the perpendicular distance of the edge from the origin. 
                    if (uniqueEdges.length==2) {
                        let APx1 = origin[0]-uniqueEdges[0][0][0];
                        let APy1 = origin[1]-uniqueEdges[0][0][1];
                        let APz1 = origin[2]-uniqueEdges[0][0][2];
                        let dx1 = uniqueEdges[0][1][0] - uniqueEdges[0][0][0];
                        let dy1 = uniqueEdges[0][1][1] - uniqueEdges[0][0][1];
                        let dz1 = uniqueEdges[0][1][2] - uniqueEdges[0][0][2];
                        let CP1 = [];
                        CP1[0] = (APy1*dz1-dy1*APz1);
                        CP1[1] = -1*(APx1*dz1-dx1*APz1);
                        CP1[2] = (APx1*dy1-dx1*APy1);
                        let D1 = Math.abs(Math.sqrt(CP1[0]**2+CP1[1]**2+CP1[2]**2)/Math.sqrt(dx1**2+dy1**2+dz1**2));
    
                        let APx2 = origin[0]-uniqueEdges[1][0][0];
                        let APy2 = origin[1]-uniqueEdges[1][0][1];
                        let APz2 = origin[2]-uniqueEdges[1][0][2];
                        let dx2 = uniqueEdges[1][1][0] - uniqueEdges[1][0][0];
                        let dy2 = uniqueEdges[1][1][1] - uniqueEdges[1][0][1];
                        let dz2 = uniqueEdges[1][1][2] - uniqueEdges[1][0][2];
                        let CP2 = [];
                        CP2[0] = (APy2*dz2-dy2*APz2);
                        CP2[1] = -1*(APx2*dz2-dx2*APz2);
                        CP2[2] = (APx2*dy2-dx2*APy2);
                        let D2 = Math.abs(Math.sqrt(CP2[0]**2+CP2[1]**2+CP2[2]**2)/Math.sqrt(dx2**2+dy2**2+dz2**2));
    
                        let closerEdge = [];
    
                        if (D1<D2) {
                            closerEdge.push(uniqueEdges[0]);
                        }
                        else {
                            closerEdge.push(uniqueEdges[1]);
                        }
                        // Appending the closer edge and the corresponding coordinates in the extreme coordinates. 
                        for (let i=0; i<intersectionCoordinateEdgesUnique.length; i++) {
                            if ((JSON.stringify(intersectionCoordinateEdgesUnique[i])==JSON.stringify(closerEdge[0])) || 
                            (JSON.stringify([intersectionCoordinateEdgesUnique[i][1],intersectionCoordinateEdgesUnique[i][0]])==JSON.stringify(closerEdge[0]))){
                                extremeCoordinateEdges.push(intersectionCoordinateEdgesUnique[i]);
                                extremeCoordinates.push(intersectionCoordinatesAtTopUnique[i]);
                            }
                        }    
                    }
                    else if (uniqueEdges.length==3){
                        // The only possible scenario in this case will be if the two points lie on two parallel edges and rest of the two lie on one edge. 
                        // The points lying on the parallel edges will always be the governing ones in this case. 
                        let sameEdgeIndex = [];
                        for (let i=0; i<intersectionCoordinatesAtTopUnique.length-1; i++) {
                            for(let j=i+1; j<intersectionCoordinatesAtTopUnique.length; j++) {
                                if (JSON.stringify(intersectionCoordinateEdgesUnique[i])==JSON.stringify(intersectionCoordinateEdgesUnique[j]) ||
                                JSON.stringify(intersectionCoordinateEdgesUnique[i])==JSON.stringify([intersectionCoordinateEdgesUnique[j][1],intersectionCoordinateEdgesUnique[j][0]])){
                                    sameEdgeIndex.push(i);
                                    sameEdgeIndex.push(j);
                                }
                            }
                        }
                        sameEdgeIndex = [...new Set(sameEdgeIndex)];
                        // The sameEdgeIndex will contain the indices from intersectionCoordinateEdgesUnique, which are repreated. 
                        // the edges other than these should be essentially the other two parallel edges. 
                        for (let i=0; i<intersectionCoordinateEdgesUnique.length; i++) {
                            if (!sameEdgeIndex.includes(i)) {
                                extremeCoordinateEdges.push(intersectionCoordinateEdgesUnique[i]);
                                extremeCoordinates.push(intersectionCoordinatesAtTopUnique[i]);
                            }
                        }
                    }
                }
                else if (intersectionCoordinatesAtTopUnique.length > 4) {
                // In this case, the extreme points will always lie on the surface closer to the origin. 
                let distance = [];
                for (let i=0; i<intersectionCoordinatesAtTopUnique.length; i++) {
                    distance.push(getLineDistanceFromPoint(intersectionCoordinateEdgesUnique[i], origin));
                }
                let minimumDistance = math.min(...distance);
                let minimumDistanceIndex = [];
                // minimumDistanceIndex will contain the indices of the surfaces which have minimum value of distance from the origin. 
                // Then we can append the coordinates and corresponding surfaces based on these indices. 
                for (let i=0; i<distance.length; i++) {
                    if (distance[i]===minimumDistance) {
                        minimumDistanceIndex.push(i);
                    }
                }
                for (let i=0; i<minimumDistanceIndex.length; i++){
                    extremeCoordinates.push(intersectionCoordinatesAtTopUnique[minimumDistanceIndex[i]]);
                    extremeCoordinateSurfaces.push(intersectionCoordinateEdgesUnique[[minimumDistanceIndex[i]]]);
                } 
                }
                else {
                    extremeCoordinateEdges = intersectionCoordinateEdgesUnique;
                    extremeCoordinates = intersectionCoordinatesAtTopUnique;
                } 
                return {extremeCoordinateEdges,extremeCoordinates,coordinateOn:"edge"};
            }
            else if (heightOfRotatingObject < heightPolyhedralSurface){
                // This is the case where the polyhedral surface is higher than the triangulated surface.
                // The key here, in my understanding, is that the extreme points will be the ones lying on the traingle edges which lie on the circumference of the circle. 
                // So, we need to find out those edges and the corresponding points of intersection. 
                // According to my understanding, there wil be similar situations as in the previous case. The possible intersection coordinates will be 4,3,2,1. 
                for (let i=0; i<coordinateOfIntersectionEdges.length; i++) {
                    // Check if the two z values are same as the height of the crane.
                    // You also need to add a check that both the points lie on the outer circle. Otherwise, you will get the coordinates from inside as well. 
                    // You do not need to check the possibility of coordinates lying on the polyhedral surface edges here because the height of polyhedral surface is higher than the crane.
                    // So here we are finding only those intersection coordinates which lie on the edges on the boundary. 
                    if ((coordinateOfIntersectionEdges[i][0][2].toFixed(4)===coordinateOfIntersectionEdges[i][1][2].toFixed(4)) 
                    && (parseFloat(coordinateOfIntersectionEdges[i][0][2].toFixed(4))===heightOfRotatingObject)
                    && (parseFloat(Math.sqrt((coordinateOfIntersectionEdges[i][0][0]-origin[0])**2+(coordinateOfIntersectionEdges[i][0][1]-origin[1])**2).toFixed(4))===radiusCircle)
                    && (parseFloat(Math.sqrt((coordinateOfIntersectionEdges[i][1][0]-origin[0])**2+(coordinateOfIntersectionEdges[i][1][1]-origin[1])**2).toFixed(4))===radiusCircle)){
                        intersectionCoordinatesAtTop.push(coordinateOfIntersection[i]);
                        intersectionCoordinateEdges.push(coordinateOfIntersectionEdges[i]);
                    }
                }
                // Need to remove duplicates from the arrays prepared above. This is because the same edge might occur multiple times as it might be common in different cube faces. 
                // Another thing is edge coordinates may alse be reverse but same. This also has to be taken care off.
                for (let i=0; i<intersectionCoordinatesAtTop.length; i++) {
                    if (i==0){
                        intersectionCoordinatesAtTopUnique.push(intersectionCoordinatesAtTop[i]);
                        intersectionCoordinateEdgesUnique.push(intersectionCoordinateEdges[i]);
                    }
                    if (i!=0){
                        let count = 0;
                        for (let j=0; j<i; j++) {
                            if ((JSON.stringify(intersectionCoordinatesAtTop[i])==JSON.stringify(intersectionCoordinatesAtTop[j])) && 
                            ((JSON.stringify(intersectionCoordinateEdges[i])==JSON.stringify(intersectionCoordinateEdges[j])) ||
                            (JSON.stringify(intersectionCoordinateEdges[i])==JSON.stringify([intersectionCoordinateEdges[j][1],intersectionCoordinateEdges[j][0]])))) {
                                count+=1;
                            }
                        }
                        if (count==0) {
                            intersectionCoordinatesAtTopUnique.push(intersectionCoordinatesAtTop[i]);
                            intersectionCoordinateEdgesUnique.push(intersectionCoordinateEdges[i]);
                        }
                    } 
                }
                // Again, similar to the previous cases, you need to find out the extreme coordinates. 
                // One way I can think off is the find out the planes of the polyhedral surface on which the points lie. 
                // One thing is, the points you might also get based on the intersection of the triangulation of the top surface of the cone and the polyhedral surface. 
                // It can give you a lot of inner coordinates easily. But, still the issue will be finding the extreme coordinates. So the better way is to use the planes I think. 
                // But this also mean that you will not be having the conditions like the intersectioncoordinates at the top are just 4 or 2 or 1. 
                // Many more coordinates will come as a result of intersection from the triangles at the top. 
                // Then check the points according to the plane.
                let intersectionCoordinatesOnSurfaces = [];
                let intersectionCoordinatesOnSurfacesEdges = [];
                let surfacesOfIntersectionCoordinatesOnSurfaces = [];
                // This loop should result in three arrays, the coordinates which lie on the polyhedral surfaces, the triangulated surface edges on which they lie and the surfaces themselves which they satisfy. 
                for (let k=0; k<intersectionCoordinatesAtTopUnique.length; k++) {
                    for (let l=0; l<polyhedralSurface.length; l++) {
                        let planeEquation = getSurfaceEquation(polyhedralSurface[l]);
                        let ifPointLiesOnPlane = checkIfPointLiesOnPlane(intersectionCoordinatesAtTopUnique[k], planeEquation);
                        if (ifPointLiesOnPlane==true) {
                            let ifPointLiesWithinLimits = checkIfLiesWithinLimit(polyhedralSurface[l], intersectionCoordinatesAtTopUnique[k], "TriangleEdgetoCubeSurface");
                            if (ifPointLiesWithinLimits==true) {
                                intersectionCoordinatesOnSurfaces.push(intersectionCoordinatesAtTopUnique[k]);
                                intersectionCoordinatesOnSurfacesEdges.push(intersectionCoordinateEdgesUnique[k]);
                                surfacesOfIntersectionCoordinatesOnSurfaces.push(polyhedralSurface[l]);
                            }
                        }
                    }
                }
                if (intersectionCoordinatesOnSurfaces.length==4) {
                    let uniqueSurfaces = [];
                    // Getting all the unique surfaces, which contain the intersection points. 
                    uniqueSurfaces.push(surfacesOfIntersectionCoordinatesOnSurfaces[0]);
                    for (let i=1; i<surfacesOfIntersectionCoordinatesOnSurfaces.length; i++){
                        let count = 0;
                        for (let j=0; j<i; j++) {
                            if (JSON.stringify(surfacesOfIntersectionCoordinatesOnSurfaces[i])===JSON.stringify(surfacesOfIntersectionCoordinatesOnSurfaces[j])) {
                                count=count+1;
                            }
                        }
                        if (count==0) {
                            uniqueSurfaces.push(surfacesOfIntersectionCoordinatesOnSurfaces[i]);
                        }
                    }
                    //If the length of the unique surfaces is equal to two, this means that 4 intersection points are lying on 2 parallel surfaces. 
                    // The surface closer to the origin should be the governing one in this case. 
                    if (uniqueSurfaces.length==2){
                        // Find the distance of the origin from the planes. 
                        let plane1 = getSurfaceEquation(uniqueSurfaces[0]);
                        let plane2 = getSurfaceEquation(uniqueSurfaces[1]);
                        let distance1 = Math.abs(plane1.normalVector0*origin[0]+plane1.normalVector1*origin[1]+plane1.normalVector2*origin[2]+plane1.d)/Math.sqrt(plane1.normalVector0**2+plane1.normalVector1**2+plane1.normalVector2**2);
                        let distance2 = Math.abs(plane2.normalVector0*origin[0]+plane2.normalVector1*origin[1]+plane2.normalVector2*origin[2]+plane2.d)/Math.sqrt(plane2.normalVector0**2+plane2.normalVector1**2+plane2.normalVector2**2);
                        if (distance1<distance2) {
                            var closerSurface = uniqueSurfaces[0];
                        }
                        else if (distance2<distance1) {
                            var closerSurface = uniqueSurfaces[1];
                        }
                        for (let i=0; i<surfacesOfIntersectionCoordinatesOnSurfaces.length; i++){
                            if (JSON.stringify(surfacesOfIntersectionCoordinatesOnSurfaces[i])===JSON.stringify(closerSurface)) {
                                extremeCoordinates.push(intersectionCoordinatesOnSurfaces[i]);
                                extremeCoordinateSurfaces.push(surfacesOfIntersectionCoordinatesOnSurfaces[i]);
                            }
                        }
                    }
                    else if (uniqueSurfaces.length==3) {
                        // One possible case in this is when the intersection occurs between two parallel surfaces and one non parallel face. 
                        // So, finding the points lying on the parallel surfaces should do the trick. 
                        let sameSurfaceIndex = [];
                        for (let i=0; i<intersectionCoordinatesOnSurfaces.length-1; i++) {
                            for(let j=i+1; j<intersectionCoordinatesOnSurfaces.length; j++) {
                                if (JSON.stringify(surfacesOfIntersectionCoordinatesOnSurfaces[i])===JSON.stringify(surfacesOfIntersectionCoordinatesOnSurfaces[j])) {
                                    sameSurfaceIndex[0] = i;
                                    sameSurfaceIndex[1] = j;
                                    break;
                                }
                            }
                        }
                        // sameSurfaceIndex will contain the indices of the points which are lying on the same surface. 
                        for (let i=0; i<intersectionCoordinatesOnSurfaces.length; i++) {
                            if (i!=sameSurfaceIndex[0] && i!=sameSurfaceIndex[1]) {
                                extremeCoordinates.push(intersectionCoordinatesOnSurfaces[i]);
                                extremeCoordinateSurfaces.push(surfacesOfIntersectionCoordinatesOnSurfaces[i]);
                            }
                        }
                        // Another scenario in this case can be if one of the intersection point lies on the corner of the polyhedral surface. 
                        // That, I will handle later
                    }
                }
                else if (intersectionCoordinatesOnSurfaces.length > 4) {
                    // In this case, the extreme points will always lie on the surface closer to the origin. 
                    let distance = [];
                    for (let i=0; i<intersectionCoordinatesOnSurfaces.length; i++) {
                        distance.push(getPlaneDistanceFromPoint(surfacesOfIntersectionCoordinatesOnSurfaces[i], origin));
                    }
                    let minimumDistance = math.min(...distance);
                    let minimumDistanceIndex = [];
                    // minimumDistanceIndex will contain the indices of the surfaces which have minimum value of distance from the origin. 
                    // Then we can append the coordinates and corresponding surfaces based on these indices. 
                    for (let i=0; i<distance.length; i++) {
                        if (distance[i]===minimumDistance) {
                            minimumDistanceIndex.push(i);
                        }
                    }
                    for (let i=0; i<minimumDistanceIndex.length; i++){
                        extremeCoordinates.push(intersectionCoordinatesOnSurfaces[minimumDistanceIndex[i]]);
                        extremeCoordinateSurfaces.push(surfacesOfIntersectionCoordinatesOnSurfaces[[minimumDistanceIndex[i]]]);
                    }
                }
                else {
                    extremeCoordinates = intersectionCoordinatesOnSurfaces;
                    extremeCoordinateSurfaces = surfacesOfIntersectionCoordinatesOnSurfaces;
                }
                return {extremeCoordinateSurfaces,extremeCoordinates,coordinateOn:"surfaces"};
            }
        }
    
        function getClosestPoint(point1, point2, origin) {
            // So here, you need to find the point on the line joining point1 and point2, which is closest to origin in the plane of the height of the intersection point. 
            // The logic here is that we first find out the dot product of the two vectors (vector1 is point1 to origin and vector2 is point1 to point2).
            // Then we find the dot product of the two vectors. 
            // If the dot product is greater than the dot product of the base line itself, then the point will be closer to the point2.
            // Similarly, if the dot product is lesser than 0, then the point will be closer to the point1. 
            // If the dot product is lesser than the dot product of the base line itself, then the point will be in between point1 and point2.
            let adjustedOrigin = [origin[0], origin[1], point1[2]];
            let line1 = getLineEquation([point1, adjustedOrigin]);
            let line2 = getLineEquation([point1, point2]);
            let dotProduct12 = line1.a*line2.a + line1.b*line2.b + line1.c*line2.c ;
            let dotProduct22 = line2.a*line2.a + line2.b*line2.b + line2.c*line2.c ;
            if (dotProduct12>=dotProduct22) {
                return point2;
            }
            else if (dotProduct12<=0) {
                return point1;
            }
            else if (dotProduct12>0 && dotProduct12<dotProduct22){
                let projectionRatio = dotProduct12/dotProduct22;
                let closestPoint = [];
                closestPoint.push(point1[0]+line2.a*projectionRatio);
                closestPoint.push(point1[1]+line2.b*projectionRatio);
                closestPoint.push(point1[2]+line2.c*projectionRatio);
                return closestPoint;
            }
        }
    
        function getInnermostCoordinates (extremeCoordinate, origin, heightOfPartOfRotation) {
            // Getting the innermost coordinates based on the identified extreme coordinates. 
            // Innermost coordinates will actually define what needs to be cleared by the crane boom.
            // So first case is when the coordinates of intersection are on the edges. 
            // Creating a new array called innermostCoordinates to store the coordinates which are to be cleared by the crane boom. 
            let innermostCoordinates = [];
            let uniqueExtremeCoordinates = [];
            let uniqueExtremeCoordinateEdges = [];
            let uniqueExtremeCoordinatesSurfaces = [];
            if (extremeCoordinate.coordinateOn == "edge") {
                if (extremeCoordinate.extremeCoordinates.length == 1) {
                    innermostCoordinates.push(extremeCoordinate.extremeCoordinates[0]);
                }
                else if (extremeCoordinate.extremeCoordinates.length >1) {
                    // You need to first remove duplicates here. Duplicates here would be the same point lying on two different edges (if the intersection point is the corner of the surfaces).
                    uniqueExtremeCoordinates.push(extremeCoordinate.extremeCoordinates[0]);
                    uniqueExtremeCoordinateEdges.push(extremeCoordinate.extremeCoordinateEdges[0]);
                    for (let i=1; i<extremeCoordinate.extremeCoordinates.length; i++){
                        let count = 0;
                        for (let j=0; j<i; j++){
                            if (JSON.stringify(extremeCoordinate.extremeCoordinates[i])==JSON.stringify(extremeCoordinate.extremeCoordinates[j])) {
                                count = count+1;
                            }
                        }
                        if (count==0) {
                            uniqueExtremeCoordinates.push(extremeCoordinate.extremeCoordinates[i]);
                            uniqueExtremeCoordinateEdges.push(extremeCoordinate.extremeCoordinateEdges[i]);
                        }
                    }
                    // First case, when the unique intersection coordinate is only one. This means that it is just touching the surface.
                    if (uniqueExtremeCoordinates.length == 1) {
                        innermostCoordinates.push(uniqueExtremeCoordinates[0]);
                    }
                    else if (uniqueExtremeCoordinates.length == 2){
                        // Check if the two lines are parallel 
                        let ifParallel = checkIfLinesAreparallel (uniqueExtremeCoordinateEdges[0], uniqueExtremeCoordinateEdges[1]);
                        // First case is when the two lines are not parallel. This means they will intersect with each other.
                        // In this case, the point of intersection will be innermost point. 
                        // You can also add the coordinates of the midpoint of the intersection point and the extreme coordinates.
                        // interestingly, the intersection coordinate in this case would be the point which is common in both the edges. You don't need to find intersection separately.  
                        if (ifParallel == false){
                            for (let i=0; i<uniqueExtremeCoordinateEdges[0].length; i++) {
                                for (let j=0; j<uniqueExtremeCoordinateEdges[1].length; j++){
                                    if (JSON.stringify(uniqueExtremeCoordinateEdges[0][i]) == JSON.stringify(uniqueExtremeCoordinateEdges[1][j])){
                                        innermostCoordinates.push(uniqueExtremeCoordinateEdges[0][i]);
                                    }
                                }
                            }
                            let average1 = [(uniqueExtremeCoordinates[0][0]+innermostCoordinates[0][0])/2, (uniqueExtremeCoordinates[0][1]+innermostCoordinates[0][1])/2, (uniqueExtremeCoordinates[0][2]+innermostCoordinates[0][2])/2];
                            let average2 = [(uniqueExtremeCoordinates[1][0]+innermostCoordinates[0][0])/2, (uniqueExtremeCoordinates[1][1]+innermostCoordinates[0][1])/2, (uniqueExtremeCoordinates[1][2]+innermostCoordinates[0][2])/2];
                            // The centre point might not always be the closest coordinates. We need to find the point on the line which is closest to the centre point. 
                            let closestPoint1 = getClosestPoint(uniqueExtremeCoordinates[0], innermostCoordinates[0], origin);
                            let closestPoint2 = getClosestPoint(uniqueExtremeCoordinates[1], innermostCoordinates[0], origin);
                            innermostCoordinates.push(average1);
                            innermostCoordinates.push(average2);
                            innermostCoordinates.push(closestPoint1);
                            innermostCoordinates.push(closestPoint2);
                        }
                        else {
                            // In case when they are parallel, there can be two cases. 
                            // One, when the points lie on one edge, i.e. both the lines coincide. 
                            // Second the points lie on two different parallel edges. 
                            if ((JSON.stringify(uniqueExtremeCoordinateEdges[0])==JSON.stringify(uniqueExtremeCoordinateEdges[1])) ||
                                (JSON.stringify(uniqueExtremeCoordinateEdges[0])==JSON.stringify([uniqueExtremeCoordinateEdges[1][1], uniqueExtremeCoordinateEdges[1][0]]))) {
                                    // This is the case where both the points lie on the same edge. 
                                    // The innermost coordinates in this case would be the ones the intersection coordinates themselves. 
                                    // Also, take the average into account. 
                                    innermostCoordinates.push([(uniqueExtremeCoordinates[0][0]+uniqueExtremeCoordinates[1][0])/2, (uniqueExtremeCoordinates[0][1]+uniqueExtremeCoordinates[1][1])/2, (uniqueExtremeCoordinates[0][2]+uniqueExtremeCoordinates[1][2])/2]);
                                    innermostCoordinates.push(uniqueExtremeCoordinates[0]);
                                    innermostCoordinates.push(uniqueExtremeCoordinates[1]);
                                    let closestPoint1 = getClosestPoint(uniqueExtremeCoordinates[0], uniqueExtremeCoordinates[1], origin);
                                    innermostCoordinates.push(closestPoint1);
                                }
                            else {
                                // This will be the case when the coordinates lie on the two different parallel edges. 
                                // First, check if both the coordinates have +ve x values. The innermost point will be left side end of the edge.
                                // Simply speaking, the edges closer to the origin will be the innermost ones.  
                                if ((uniqueExtremeCoordinates[0][0] > 0) && (uniqueExtremeCoordinates[1][0] > 0)) {
                                    let closerCoordinates = getCloserCoordinates(uniqueExtremeCoordinateEdges, origin);
                                    let averageOfCloserCoordinates = [(closerCoordinates[0][0]+closerCoordinates[1][0])/2, (closerCoordinates[0][1]+closerCoordinates[1][1])/2, (closerCoordinates[0][2]+closerCoordinates[1][2])/2]; 
                                    innermostCoordinates.push(closerCoordinates[0]);
                                    innermostCoordinates.push(closerCoordinates[1]);
                                    innermostCoordinates.push(averageOfCloserCoordinates);
                                    let closestPoint1 = getClosestPoint(closerCoordinates[0], closerCoordinates[1], origin);
                                    let closestPoint2 = getClosestPoint(closerCoordinates[0], uniqueExtremeCoordinates[0], origin);
                                    let closestPoint3 = getClosestPoint(closerCoordinates[1], uniqueExtremeCoordinates[1], origin);
                                    innermostCoordinates.push(closestPoint1);
                                    innermostCoordinates.push(closestPoint2);
                                    innermostCoordinates.push(closestPoint3);
                                }
                                if ((uniqueExtremeCoordinates[0][0] < 0) && (uniqueExtremeCoordinates[1][0] < 0)) {
                                    let closerCoordinates = getCloserCoordinates(uniqueExtremeCoordinateEdges, origin);
                                    let averageOfCloserCoordinates = [(closerCoordinates[0][0]+closerCoordinates[1][0])/2, (closerCoordinates[0][1]+closerCoordinates[1][1])/2, (closerCoordinates[0][2]+closerCoordinates[1][2])/2]; 
                                    innermostCoordinates.push(closerCoordinates[0]);
                                    innermostCoordinates.push(closerCoordinates[1]);
                                    innermostCoordinates.push(averageOfCloserCoordinates);
                                    let closestPoint1 = getClosestPoint(closerCoordinates[0], closerCoordinates[1], origin);
                                    let closestPoint2 = getClosestPoint(closerCoordinates[0], uniqueExtremeCoordinates[0], origin);
                                    let closestPoint3 = getClosestPoint(closerCoordinates[1], uniqueExtremeCoordinates[1], origin);
                                    innermostCoordinates.push(closestPoint1);
                                    innermostCoordinates.push(closestPoint2);
                                    innermostCoordinates.push(closestPoint3);
                                }
                                // Now, check if the sign of the X is different for both the coordinates.
                                if (((uniqueExtremeCoordinates[0][0] < 0) && (uniqueExtremeCoordinates[1][0] > 0)) ||
                                ((uniqueExtremeCoordinates[0][0] > 0) && (uniqueExtremeCoordinates[1][0] < 0))) {
                                    let closerCoordinates = getCloserCoordinates(uniqueExtremeCoordinateEdges, origin);
                                    let averageOfCloserCoordinates = [(closerCoordinates[0][0]+closerCoordinates[1][0])/2, (closerCoordinates[0][1]+closerCoordinates[1][1])/2, (closerCoordinates[0][2]+closerCoordinates[1][2])/2]; 
                                    innermostCoordinates.push(closerCoordinates[0]);
                                    innermostCoordinates.push(closerCoordinates[1]);
                                    innermostCoordinates.push(averageOfCloserCoordinates);
                                    let closestPoint1 = getClosestPoint(closerCoordinates[0], closerCoordinates[1], origin);
                                    let closestPoint2 = getClosestPoint(closerCoordinates[0], uniqueExtremeCoordinates[0], origin);
                                    let closestPoint3 = getClosestPoint(closerCoordinates[1], uniqueExtremeCoordinates[1], origin);
                                    innermostCoordinates.push(closestPoint1);
                                    innermostCoordinates.push(closestPoint2);
                                    innermostCoordinates.push(closestPoint3);
                                }
                            }
                        }
                    }
                }
            }
            // The second case is when the coordinates of intersection are on the surfaces of the pulyhedral surface. This means that the crane boom height is smaller than object. 
            else if (extremeCoordinate.coordinateOn == "surfaces") {
                if (extremeCoordinate.extremeCoordinates.length == 1) {
                    innermostCoordinates.push(extremeCoordinate.extremeCoordinates[0]);
                }
                else if (extremeCoordinate.extremeCoordinates.length >1) {
                    // You need to first remove duplicates here. Duplicates here would be the same point lying on two different surfaces (if the intersection point is the corner of the surface).
                    uniqueExtremeCoordinates.push(extremeCoordinate.extremeCoordinates[0]);
                    uniqueExtremeCoordinatesSurfaces.push(extremeCoordinate.extremeCoordinateSurfaces[0]);
                    for (let i=1; i<extremeCoordinate.extremeCoordinates.length; i++){
                        let count = 0;
                        for (let j=0; j<i; j++){
                            if (JSON.stringify(extremeCoordinate.extremeCoordinates[i])===JSON.stringify(extremeCoordinate.extremeCoordinates[j])) {
                                count = count+1;
                            }
                        }
                        if (count==0) {
                            uniqueExtremeCoordinates.push(extremeCoordinate.extremeCoordinates[i]);
                            uniqueExtremeCoordinatesSurfaces.push(extremeCoordinate.extremeCoordinateSurfaces[i]);
                        }
                    }
                    // First case, when the unique intersection coordinate is only one. This means that it is just touching the surface.
                    if (uniqueExtremeCoordinates.length == 1) {
                        innermostCoordinates.push(uniqueExtremeCoordinates[0]);
                    }
                    else if (uniqueExtremeCoordinates.length == 2){
                        // Check if the two surfaces are parallel 
                        let ifParallel = checkIfSurfacesAreparallel (uniqueExtremeCoordinatesSurfaces[0], uniqueExtremeCoordinatesSurfaces[1]);
                        // First case is when the two surfaces are not parallel. This means they will intersect with each other.
                        // In this case, the point of intersection will be innermost point. 
                        // You can also add the coordinates of the midpoint of the intersection point and the extreme coordinates.
                        // interestingly, the intersection coordinate in this case would be the point which is common in both the edges. You don't need to find intersection separately.
                        // So here, you need to find the X and Y coordinates of the intersection of the surfaces.
                        // The Z coordinate will be the height of the crane/ triangulated surface.   
                        if (ifParallel == false){
                            // Get the top lines of the surfaces. As the surfaces are always perpendicular to X-Y plane, I am trying to get the lines in X-Y plane which can represent them.  
                            let height = [];
                            for (let i=0; i<uniqueExtremeCoordinatesSurfaces.length; i++) {
                                for (let j=0; j<uniqueExtremeCoordinatesSurfaces[i].length; j++) {
                                    height.push(uniqueExtremeCoordinatesSurfaces[i][j][2]);
                                }
                            }
                            let heightMax = Math.max(...height);
                            // Getting the coordinates which lie at the top. 
                            let coordinatesAtTop = [];
                            for (let i=0; i<uniqueExtremeCoordinatesSurfaces.length; i++) {
                                for (let j=0; j<uniqueExtremeCoordinatesSurfaces[i].length; j++) {
                                    if(uniqueExtremeCoordinatesSurfaces[i][j][2]==heightMax){
                                        coordinatesAtTop.push(uniqueExtremeCoordinatesSurfaces[i][j]);
                                    }
                                }
                            }
                            // Getting the distance of the coordinates at the top based on their X,y coordinates. 
                            let distanceFromOrigin = [];
                            for (let i=0; i<coordinatesAtTop.length; i++) {
                                distanceFromOrigin.push(Math.sqrt((coordinatesAtTop[i][0]-origin[0])**2+(coordinatesAtTop[i][1]-origin[1])**2));
                            }
                            let minDistanceFromOrigin = Math.min(...distanceFromOrigin);
                            // Appending the coordinates at the top which have the minimum distance from the origin. 
                            for (let i=0; i<distanceFromOrigin.length; i++){
                                if (distanceFromOrigin[i]==minDistanceFromOrigin) {
                                    innermostCoordinates.push(coordinatesAtTop[i]);
                                }
                            }
                            // Changing the z value to the value of the height of the triangulated surface. 
                            for (let i=0; i<innermostCoordinates.length; i++) {
                                innermostCoordinates[i][2] = heightOfPartOfRotation;
                            }
                            let average1 = [(uniqueExtremeCoordinates[0][0]+innermostCoordinates[0][0])/2, (uniqueExtremeCoordinates[0][1]+innermostCoordinates[0][1])/2, (uniqueExtremeCoordinates[0][2]+innermostCoordinates[0][2])/2];
                            let average2 = [(uniqueExtremeCoordinates[1][0]+innermostCoordinates[0][0])/2, (uniqueExtremeCoordinates[1][1]+innermostCoordinates[0][1])/2, (uniqueExtremeCoordinates[1][2]+innermostCoordinates[0][2])/2];
                            let closestPoint1 = getClosestPoint(uniqueExtremeCoordinates[0], innermostCoordinates[0], origin);
                            let closestPoint2 = getClosestPoint(uniqueExtremeCoordinates[1], innermostCoordinates[0], origin);
                            innermostCoordinates.push(average1);
                            innermostCoordinates.push(average2);
                            innermostCoordinates.push(closestPoint1);
                            innermostCoordinates.push(closestPoint2);
                        }
                        else {
                            // In case when they are parallel, there can be two cases. 
                            // One, when the points lie on one surface, i.e. both the surfaces coincide. 
                            // Second the points lie on two different parallel surfaces.
                            // To check if the surfaces are the same.  
                            let height = [];
                            for (let i=0; i<uniqueExtremeCoordinatesSurfaces.length; i++) {
                                for (let j=0; j<uniqueExtremeCoordinatesSurfaces[i].length; j++) {
                                    height.push(uniqueExtremeCoordinatesSurfaces[i][j][2]);
                                }
                            }
                            let heightMax = Math.max(...height);
                            // Getting the coordinates which lie at the top for each polyhedral surface. 
                            let coordinatesAtTop = [];
                            for (let i=0; i<uniqueExtremeCoordinatesSurfaces.length; i++) {
                                coordinatesAtTop[i] = [];
                                for (let j=0; j<uniqueExtremeCoordinatesSurfaces[i].length; j++) {
                                    if(uniqueExtremeCoordinatesSurfaces[i][j][2]==heightMax){
                                        coordinatesAtTop[i].push([uniqueExtremeCoordinatesSurfaces[i][j][0],uniqueExtremeCoordinatesSurfaces[i][j][1]]);
                                    }
                                }
                            }
                            let xofCoordinatesAtTop = [];
                            let yofCoordinatesAtTop = [];
                            // Getting the x and y coordinates separately in a list
                            for (let i=0; i<coordinatesAtTop.length; i++){
                                xofCoordinatesAtTop[i] = [];
                                yofCoordinatesAtTop[i] = [];
                                for (let j=0; j<coordinatesAtTop[i].length; j++) {
                                    xofCoordinatesAtTop[i].push (coordinatesAtTop[i][j][0]);
                                    yofCoordinatesAtTop[i].push (coordinatesAtTop[i][j][1]);
                                }
                            }
                            // Getting the minimum and maximum values of x and y for the coordinates at the top
                            let xminmax = [];
                            let yminmax = [];
                            for (let i=0; i<xofCoordinatesAtTop.length; i++){
                                xminmax[i] = [];
                                yminmax[i] = [];
                                xminmax[i] = [Math.min(...xofCoordinatesAtTop[i]), Math.max(...xofCoordinatesAtTop[i])];
                                yminmax[i] = [Math.min(...yofCoordinatesAtTop[i]), Math.max(...yofCoordinatesAtTop[i])];
                            }
    
                            // If minimum and maximum value of x and y are same for both the surfaces, this means the surfaces should be equal. 
                            // They are equal means that the extreme points will be the innermost points itself
                            // Also, we can add the average of the extreme points
                            if ((xminmax[0][0]==xminmax[1][0]) && (xminmax[0][1]==xminmax[1][1]) && (yminmax[0][0]==yminmax[1][0]) && (yminmax[0][1]==yminmax[1][1])) {
                                innermostCoordinates.push(uniqueExtremeCoordinates[0]);
                                innermostCoordinates.push(uniqueExtremeCoordinates[1]);
                                innermostCoordinates.push([(uniqueExtremeCoordinates[0][0]+uniqueExtremeCoordinates[1][0])/2, (uniqueExtremeCoordinates[0][1]+uniqueExtremeCoordinates[1][1])/2, (uniqueExtremeCoordinates[0][2]+uniqueExtremeCoordinates[1][2])/2]);
                            }
                            else {
                                // This will be the case when the coordinates lie on the two different parallel surfaces. 
                                // Simply speaking, the ends closer to the origin will be the innermost ones.
                                let distance = [];
                                for (let i=0; i<coordinatesAtTop.length; i++){
                                    let distance = [];
                                    for (let j=0; j<coordinatesAtTop[i].length; j++) {
                                        distance.push(Math.sqrt((coordinatesAtTop[i][j][0]-origin[0])**2 + (coordinatesAtTop[i][j][1]-origin[1])**2));
                                    }
                                    let minDistance = Math.min(...distance);
                                    let maxDistanceCoordinates = [];
                                    for (let k=0; k<distance.length; k++){
                                        if (distance[k]!=minDistance) {
                                            maxDistanceCoordinates.push([coordinatesAtTop[i][k][0],coordinatesAtTop[i][k][1],heightOfPartOfRotation]);                                    
                                        }
                                    }
                                    for (let k=0; k<distance.length; k++){
                                        if (distance[k]==minDistance) {
                                            innermostCoordinates.push([coordinatesAtTop[i][k][0],coordinatesAtTop[i][k][1],heightOfPartOfRotation]); 
                                            for (let m=0; m<maxDistanceCoordinates.length; m++) {
                                                innermostCoordinates.push(getClosestPoint([coordinatesAtTop[i][k][0],coordinatesAtTop[i][k][1],heightOfPartOfRotation], maxDistanceCoordinates[m], origin));
                                            }                                   
                                        }
                                    }
                                    // In the case when the point is not the nearest one, that is the farher point.
                                    // So, we need to find the point on the line connecting the farther and closest point, which is closer to the origin.  
                                }
                                innermostCoordinates = Array.from(new Set(innermostCoordinates.map(JSON.stringify)), JSON.parse);
                                // Adding the average coordinate as well. 
                                innermostCoordinates.push([(innermostCoordinates[0][0]+innermostCoordinates[1][0])/2,(innermostCoordinates[0][1]+innermostCoordinates[1][1])/2, (innermostCoordinates[0][2]+innermostCoordinates[1][2])/2]);
                                let closestPoint1 = getClosestPoint(innermostCoordinates[0], innermostCoordinates[1], origin);
                                innermostCoordinates.push(closestPoint1);
                            }
                        }
                    }
                }
            }
            return innermostCoordinates;
        }
    
        function getAngleFromHorizontal (point1, point2) {
            let angle = Math.atan2((point2[1] - point1[1]), (point2[0] - point1[0]));
            // This gives the angle in form of a number, which basically is in radian value.
            // Angle >0 will be in first and second quadrant. So no need to make any changes. 
            if (angle>=0) {
                return angle;
            }
            // If angle is less than 0, it means the angle is in third and fourth quadrant in clockwise direction from X axis. 
            // So, subtracting it from 2*pie will give the anticlockwise angle.
            // Using the plus sign here because the angle is already negative. 
            else {
                return (2*Math.PI + angle);
            }
        }
    
        function getAdjustedGoveringCoordinate(coordinate, origin, boomClearance){
            let line = getLineEquation([origin, coordinate]);
            let perpendicularLine_a = -1*line.b;
            let perpendicularLine_b = line.a;
            // Normalizing the normal vector
            perpendicularLine_a = perpendicularLine_a/(Math.sqrt(perpendicularLine_a * perpendicularLine_a + perpendicularLine_b * perpendicularLine_b));
            perpendicularLine_b = perpendicularLine_b/(Math.sqrt(perpendicularLine_a * perpendicularLine_a + perpendicularLine_b * perpendicularLine_b));
            // Displacement in X and Y direction are identified by multiplying the normalized perpendicular vector and boom clearance. 
            let disp_x = perpendicularLine_a * boomClearance;
            let disp_y = perpendicularLine_b * boomClearance;
            return ([coordinate[0] + disp_x, coordinate[1] + disp_y]);
        }
    
        function getTriangleAreaFromVertices (vertices, origin){
            // This function gets the area of a triangle using its vertices.
            let area = 0;
            for (let i=0; i<=vertices.length-1; i++) {
                let x1 = vertices[i][0];
                let y1 = vertices[i][1];
                let x2 = vertices[i+1][0];
                let y2 = vertices[i+1][1];
                let x3 = origin[0];
                let y3 = origin[1];
                area += 0.5 * Math.abs((x1*(y2-y3) + x2*(y3-y1) + x3*(y1-y2)));
            }
            return area;
        }
    
        function getedgeRepresentationOfArcs(radiusAngleRange, origin){
            // This function returns the edges which represent each sector of the circular geometry created after detecting the intersections. 
            let edgeNumber = 100;
            vertices = {};
            for (let i=4; i<edgeNumber; i++) {
                // For each edge number, you first get the coordinates of the ends of the edges. 
                vertices[JSON.stringify(i)] = [];
                let angleIncrement = 0;
                if ((Math.max([radiusAngleRange[0], radiusAngleRange[1]]) >= Math.PI) && ((Math.max([radiusAngleRange[0], radiusAngleRange[1]]) <= 2*Math.PI))
                && (Math.min([radiusAngleRange[0], radiusAngleRange[1]]) >=0) && (Math.min([radiusAngleRange[0], radiusAngleRange[1]]) < Math.PI)) {
                    angleIncrement = (2*Math.PI - (Math.max([radiusAngleRange[0], radiusAngleRange[1]]) - Math.min([radiusAngleRange[0], radiusAngleRange[1]])))/i ; 
                }
                else {
                    angleIncrement = (Math.max([radiusAngleRange[0], radiusAngleRange[1]]) - Math.min([radiusAngleRange[0], radiusAngleRange[1]]))/i ;
                }
                for (let j=0; j<i; j++) {
                    if ((Math.max([radiusAngleRange[0], radiusAngleRange[1]]) >= Math.PI) && ((Math.max([radiusAngleRange[0], radiusAngleRange[1]]) <= 2*Math.PI))
                    && (Math.min([radiusAngleRange[0], radiusAngleRange[1]]) >=0) && (Math.min([radiusAngleRange[0], radiusAngleRange[1]]) < Math.PI)) {
                        let angle = Math.max([radiusAngleRange[0], radiusAngleRange[1]]) + j*angleIncrement;
                        let vertexX = centerX + radius * Math.cos(angle);
                        let vertexY = centerY + radius * Math.sin(angle);
                        vertices[JSON.stringify(i)].push([vertexX, vertexY]);
                    }
                    else {
                        let angle = Math.min([radiusAngleRange[0], radiusAngleRange[1]]) + j*angleIncrement;
                        let vertexX = centerX + radius * Math.cos(angle);
                        let vertexY = centerY + radius * Math.sin(angle);
                        vertices[JSON.stringify(i)].push([vertexX, vertexY]);
                    }
                }
            }
            // Finding the original sector area
            let angleDifference = 0;
            if ((Math.max([radiusAngleRange[0], radiusAngleRange[1]]) >= Math.PI) && ((Math.max([radiusAngleRange[0], radiusAngleRange[1]]) <= 2*Math.PI))
                && (Math.min([radiusAngleRange[0], radiusAngleRange[1]]) >=0) && (Math.min([radiusAngleRange[0], radiusAngleRange[1]]) < Math.PI)) {
                    angleDifference = 2*Math.PI - (Math.max([radiusAngleRange[0], radiusAngleRange[1]]) - Math.min([radiusAngleRange[0], radiusAngleRange[1]])); 
                }
            else {
                angleDifference = Math.max([radiusAngleRange[0], radiusAngleRange[1]]) - Math.min([radiusAngleRange[0], radiusAngleRange[1]]);
                }
            let originalArea = Math.PI * radiusAngleRange[2]**2 * angleDifference/(2*Math.PI);
            // Finding the area using the identified edges
            let areaEdges = [];
            let areaValues = [];
            for (let key in vertices) {
                    areaEdges.push(key);
                    areaValues.push(getTriangleAreaFromVertices (vertices[key], origin));
            }
            let areaDifference = [];
            // This list will store the values of difference in the area between the area using sector and area using edges.
            for (let i=0; i<areaValues.length; i++) {
                areaDifference.push(areaValues[i] - originalArea);
            }
            // Now the the point is to find the area where the difference between the two is not too much. 
            // This threshold will be controvertial, definitely. 
            let requiredAreaIndex = 0;
            for (let k=0; k<areaDifference.length-1; k++) {
                if ((areaDifference[k]-areaDifference[k+1])<500) {
                    requiredAreaIndex = k+1 ;
                    break;
                }
            }
            // Getting the corresponding number of edges
            let requiredEdgeNumber = 0;
            requiredEdgeNumber = areaEdges[requiredAreaIndex];
    
            let requiredVertices = vertices[requiredEdgeNumber];
            let edges = []
            for (let i=0; i<requiredVertices.length-1; i++) {
                edges.push([requiredVertices[i], requiredVertices[i+1]]);
            }
            // Appending the edges connecting the extreme coordinates and the origin. 
            // edges.push([requiredVertices[0], [origin[0],origin[1]]]);
            // edges.push([requiredVertices[requiredVertices.length-1], [origin[0],origin[1]]]);
    
            return edges; 
        }
    
        function getTangentEquation(pointOnCircle, pointCentre) {
            let slopeLine = (pointOnCircle[1] - pointCentre[1]) / (pointOnCircle[0] - pointCentre[0]);
            let tangentSlope = -1 / slopeLine;
            let constantTerm = pointOnCircle[1] - tangentSlope * pointOnCircle[0];
            return {
                slope: tangentSlope,
                intercept: constantTerm,
                point: pointOnCircle,
            }
        }
    
        function solveLineEquations(line1, line2) {
            let slope1 = line1.slope;
            let intercept1 = line1.intercept;
            let slope2 = line2.slope;
            let intercept2 = line2.intercept;
            // Calculate the x-coordinate of the intersection point
            let x = 0;
            let y = 0;
            if (!isFinite(slope2)) {
                x = line2.point[0];
                y = slope1 * x + intercept1;
            }
            else {
                x = (intercept2 - intercept1) / (slope1 - slope2);
                // Calculate the y-coordinate of the intersection point
                y = slope1 * x + intercept1;
                // Return the coordinates of the intersection point
            }
            
            return [x, y];
        }
    
        function getThetaExtents (extremeCoordinates, innermostCoordinates, radius, origin, boomClearance) {
            // This function gets the extent of angles from the horizontal line from the centre. 
            // The extent basically defines the sector which will be cut because of the presence of the obstacle.
            // The logic here is for all the points identified in the extreme and innermost coordinates, the points which have the maximum and minimum angle of rotation from the origin in anticlockwise direction
            // will be the two governing coordinates. 
            var extentGoverningCoordinates = [];
            var combinedGoverningCoordinates = extremeCoordinates.extremeCoordinates.concat(innermostCoordinates);
            var height = parseFloat(combinedGoverningCoordinates[0][2].toFixed(4));
            var adjustedOrigin = [origin[0], origin[1], height];
            // var horizontalLine = getLineEquation(adjustedOrigin, [origin[0]+1, origin[1], height]);
            var angle = [];
            for(let i=0; i<combinedGoverningCoordinates.length; i++) {
                var line = getLineEquation([adjustedOrigin, combinedGoverningCoordinates[i]]);
                angle.push(getAngleFromHorizontal(adjustedOrigin, combinedGoverningCoordinates[i]));
            }
                // Getting the maximum and minimum value indices
            let maxIndex = angle.indexOf(Math.max(...angle));
            let minIndex = angle.indexOf(Math.min(...angle));
            
            if ((angle[minIndex]>=0 && angle[minIndex] <= Math.PI) && (angle[maxIndex]>=1.5*Math.PI && angle[maxIndex] <= 2*Math.PI)) {
                // This is a slightly different case. 
                // In this case, the governing coordinates will be the ones which have the maximum angle in the first quadrant and minimum angle in the fourth quadrant. 
                let firstQuadrantAngles = [];
                let fourthQuadrantAngles = [];
                for (let j=0; j<combinedGoverningCoordinates.length; j++) {
                    if (combinedGoverningCoordinates[j][0]>=origin[0] && combinedGoverningCoordinates[j][1]>=origin[1]) {
                        // This means that the point is in first quadrant
                        firstQuadrantAngles.push (getAngleFromHorizontal(adjustedOrigin, combinedGoverningCoordinates[j]));
                    }
                    else {
                        // This means that the point is in fourth quadrant
                        fourthQuadrantAngles.push (getAngleFromHorizontal(adjustedOrigin, combinedGoverningCoordinates[j]));
                    }
                }
                maxIndex = angle.indexOf(Math.min(...fourthQuadrantAngles));
                minIndex = angle.indexOf(Math.max(...firstQuadrantAngles));
                extentGoverningCoordinates.push(combinedGoverningCoordinates[maxIndex]);
                extentGoverningCoordinates.push(combinedGoverningCoordinates[minIndex]);
            }
            else {
                extentGoverningCoordinates.push(combinedGoverningCoordinates[maxIndex]);
                extentGoverningCoordinates.push(combinedGoverningCoordinates[minIndex]);
            }
            // Getting the adjusted coordinates by adding the boom clearance.
            let distanceFromOrigin = []; 
            for (let i=0; i<extentGoverningCoordinates.length; i++) {
                // We get the distance of the extentGoverningPoints.
                distanceFromOrigin.push(Math.sqrt((extentGoverningCoordinates[i][0]-origin[0])**2 + (extentGoverningCoordinates[i][1]-origin[1])**2));
            }
            // The one which is closest to the origin should govern the adjustment of the theta.
            let minimumIndex = distanceFromOrigin.indexOf(Math.min(...distanceFromOrigin));
            let extentGoverningAdjustedCoordinate = getAdjustedGoveringCoordinate(extentGoverningCoordinates[minimumIndex], adjustedOrigin, boomClearance);
            // Finding the angle from the horizontal of the adjusted coordinate. 
            let angleFromHorizontalAdjustedCoordinate = getAngleFromHorizontal (adjustedOrigin, extentGoverningAdjustedCoordinate);
            let angleFromHorizontalGoverningCoordinate = getAngleFromHorizontal (adjustedOrigin, extentGoverningCoordinates[minimumIndex]);
            let angleDifference = Math.abs(angleFromHorizontalAdjustedCoordinate - angleFromHorizontalGoverningCoordinate);
        
            // Getting the theta of the two coordinates. 
            // These theta values will provide the extent of adjustment
            let extentTheta = [];
            for (let i=0; i<extentGoverningCoordinates.length; i++) {
                extentTheta[i] = getAngleFromHorizontal(adjustedOrigin, extentGoverningCoordinates[i]);
            }
            // Now, we can add this angleDifference in the existing theta extents.
            let minAngle = Math.min(...extentTheta);
            let maxAngle = Math.max(...extentTheta);
            if ((maxAngle>=1.5*Math.PI && maxAngle<2*Math.PI) && 
            (minAngle>=0 && minAngle<Math.PI)){
                minAngle = minAngle + angleDifference;
                maxAngle = maxAngle - angleDifference;
            }
            else {
                minAngle = minAngle - angleDifference;
                maxAngle = maxAngle + angleDifference;
            }
            return [minAngle, maxAngle];
        }
    
        function getThetaExtentsCounterweight (surface, radius, origin, boomClearance) {
            // This function gets the extent of angles from the horizontal line from the centre. 
            // The extent basically defines the sector which will be cut because of the presence of the obstacle.
            // The logic here is for all the points identified in the extreme and innermost coordinates, the points which have the maximum and minimum angle of rotation from the origin in anticlockwise direction
            // will be the two governing coordinates. 
            // This function gets the extent of angles from the horizontal line from the centre. 
            // The extent basically defines the sector which will be cut because of the presence of the obstacle.
            // The logic here is for all the points identified in the extreme and innermost coordinates, the points which have the maximum and minimum angle of rotation from the origin in anticlockwise direction
            // will be the two governing coordinates. 
            var extentGoverningCoordinates = [];
            var combinedGoverningCoordinates = surface;
            var height = parseFloat(combinedGoverningCoordinates[0].toFixed(4));
            var adjustedOrigin = [origin[0], origin[1], height];
            // var horizontalLine = getLineEquation(adjustedOrigin, [origin[0]+1, origin[1], height]);
            var angle = [];
            for(let i=0; i<combinedGoverningCoordinates.length; i++) {
                var line = getLineEquation([adjustedOrigin, combinedGoverningCoordinates[i]]);
                angle.push(getAngleFromHorizontal([adjustedOrigin, combinedGoverningCoordinates[i]]));
            }
                // Getting the maximum and minimum value indices
            let maxIndex = angle.indexOf(Math.max(...angle));
            let minIndex = angle.indexOf(Math.min(...angle));
            
            if ((minIndex>=0 && minIndex <= Math.PI) && (maxIndex>=1.5*Math.PI && MaxIndex <= 2*Math.PI)) {
                // This is a slightly different case. 
                // In this case, the governing coordinates will be the ones which have the maximum angle in the first quadrant and minimum angle in the fourth quadrant. 
                let firstQuadrantAngles = [];
                let fourthQuadrantAngles = [];
                for (let j=0; j<combinedGoverningCoordinates.length; j++) {
                    if (combinedGoverningCoordinates[j][0]>=origin[0] && combinedGoverningCoordinates[j][1]>=origin[1]) {
                        // This means that the point is in first quadrant
                        firstQuadrantAngles.push (getAngleFromHorizontal([adjustedOrigin, combinedGoverningCoordinates[j]]));
                    }
                    else {
                        // This means that the point is in fourth quadrant
                        fourthQuadrantAngles.push (getAngleFromHorizontal([adjustedOrigin, combinedGoverningCoordinates[j]]));
                    }
                }
                maxIndex = angle.indexOf(Math.min(...fourthQuadrantAngles));
                minIndex = angle.indexOf(Math.max(...firstQuadrantAngles));
                extentGoverningCoordinates.push(combinedGoverningCoordinates[maxIndex]);
                extentGoverningCoordinates.push(combinedGoverningCoordinates[minIndex]);
            }
            else {
                extentGoverningCoordinates.push(combinedGoverningCoordinates[maxIndex]);
                extentGoverningCoordinates.push(combinedGoverningCoordinates[minIndex]);
            }
            // Getting the adjusted coordinates by adding the boom clearance. 
            for (let i=0; i<extentGoverningCoordinates.length; i++) {
                extentGoverningCoordinates[i] = getAdjustedGoveringCoordinate(extentGoverningCoordinates[i], adjustedOrigin, boomClearance);
            }
            // Getting the theta of the two coordinates. 
            // These theta values will provide the extent of adjustment
            var extentTheta = [];
            for (let i=0; i<extentGoverningCoordinates.length; i++) {
                extentTheta[i] = getAngleFromHorizontal(adjustedOrigin, extentGoverningCoordinates[i]);
            }
            return extentTheta;
        }
    
        function getAdjustedRadius (boomLength, innermostCoordinates, origin, boomClearance, heightOfPolyhedralSurface, heightOfObjectOfRotation) {
            // Finding the point from the innermost coordinates which is closest to the origin. 
            // Need to take the boom clearance into consideration as well. 
            let dist = [];
            for (let i=0; i<innermostCoordinates.length; i++) {
                dist.push(Math.sqrt((innermostCoordinates[i][0]-origin[0])**2 + (innermostCoordinates[i][1]-origin[1])**2));
            }
            let minValue = Math.min(...dist);
            let minIndex = dist.indexOf(minValue);
            let governingCoordinate = innermostCoordinates[minIndex];
            if (heightOfPolyhedralSurface<=heightOfObjectOfRotation) {
                let Dt0 = minValue;
                let H0 = governingCoordinate[2];
                let theta1 = Math.atan(H0/Dt0);
                let theta2 = Math.asin(boomClearance/Math.sqrt(Dt0**2 + H0**2));
                let deltaH = boomClearance/Math.cos(theta1 + theta2);
                let Radj = Dt0*boomLength / (Math.sqrt(Dt0**2 + (H0+deltaH)**2));
                return Radj;
            }
            if (heightOfPolyhedralSurface>heightOfObjectOfRotation) {
                // In this case, there will be two situations.
                // First, when after the adjustment, the crane boom passes over the polyhedral surface
                let Dt0 = minValue;
                let H0 = heightOfPolyhedralSurface;
                if (Math.sqrt(Dt0**2 + H0**2) <= boomLength) {
                    let theta1 = Math.atan(H0/Dt0);
                    let theta2 = Math.asin(boomClearance/Math.sqrt(Dt0**2 + H0**2));
                    let deltaH = boomClearance/Math.cos(theta1 + theta2);
                    let Radj = Dt0*boomLength / (Math.sqrt(Dt0**2 + (H0+deltaH)**2));
                    return Radj;
                }
                else if(Math.sqrt(Dt0**2 + H0**2) > boomLength) {
                    // Second, when the crane boom length is not enough to pass over the polyhedral surface.
                    // In this case, the crane has to clear the height which formed based on the angle made by maximum boom length.
                    // Simpley speaking, in this case, Radj has to be the Distance of the object from crane centre minus the boom clearance.  
                    //let theta = Math.cos()
                    //let theta1 = Math.atan(H0/Dt0);
                    //let theta2 = Math.asin(boomClearance/Math.sqrt(Dt0**2 + H0**2));
                    //let deltaH = boomClearance/Math.cos(theta1 + theta2);
                    let Radj = Dt0 - boomClearance;
                    return Radj;
                }   
            } 
        } 
    
        function getAdjustedGeometryPolar (thetaExtent, adjustedRadius, radius, origin) {
            // Sorting the theta extents based on the first value of of each inner list.
            let extentTheta =  thetaExtent;
            let sortedThetaExtent = extentTheta.sort((a, b) => a[0] - b[0]);
            let radiusAngleRanges = [];
            // So here, for all the values in the sortedThetaExtent, you store the radius along with it. 
            // First storing for those values which have the radius less than the outer radius.
            for (let i=0; i<sortedThetaExtent.length; i++) {
                for (let j=0; j<thetaExtent.length; j++) {
                    if (JSON.stringify(sortedThetaExtent[i]) === JSON.stringify(thetaExtent[j])) {
                        radiusAngleRanges.push([sortedThetaExtent[i][0], sortedThetaExtent[i][1], adjustedRadius[j]]);
                    }
                }
            }
            // Then storing the values for which the radius is equal to the outer radius. 
            for (let i=0; i<sortedThetaExtent.length-1; i++) {
                radiusAngleRanges.push([sortedThetaExtent[i][1], sortedThetaExtent[i+1][0], radius]);
            }
            // Adding the value corresponding to the last sortedThetaExtent and the first one. 
            radiusAngleRanges.push([sortedThetaExtent[sortedThetaExtent.length-1][1], sortedThetaExtent[0][0], radius]);
            return radiusAngleRanges;
        }
    
        // The next thing we need to see is if the lifted object is colliding with something. 
        // The cylinder of rotation of the lifted object should be developed at the set elevation. 
        // Then you check the intersection between the cylinder and the objects. 
        // Based on the intersection, you need to prepare the updated geometry. 
        function getCSpace (surfaceAtMaximumHeight, liftedObjectRadius) {
            // Here, you need to create the Configuration Space (C-Space) for the surface. 
            let newCoordinates = [];
            for (let i=0; i<surfaceAtMaximumHeight.length-1; i++) {
                let newPoints = [];
                // Calculate the direction vector of the line segments joining the current point to its previous and next point.
                if (i!=0) {
                    var direction1 = [surfaceAtMaximumHeight[i+1][0] - surfaceAtMaximumHeight[i][0], surfaceAtMaximumHeight[i+1][1] - surfaceAtMaximumHeight[i][1]];
                    var direction2 = [surfaceAtMaximumHeight[i][0] - surfaceAtMaximumHeight[i-1][0], surfaceAtMaximumHeight[i][1] - surfaceAtMaximumHeight[i-1][1]];
                }
                else {
                    var direction1 = [surfaceAtMaximumHeight[i+1][0] - surfaceAtMaximumHeight[i][0], surfaceAtMaximumHeight[i+1][1] - surfaceAtMaximumHeight[i][1]];
                    var direction2 = [surfaceAtMaximumHeight[i][0] - surfaceAtMaximumHeight[surfaceAtMaximumHeight.length-2][0], surfaceAtMaximumHeight[i][1] - surfaceAtMaximumHeight[surfaceAtMaximumHeight.length-2][1]];
                }
                
                // Calculate the length of the line segment
                let length1 = Math.sqrt(direction1[0] ** 2 + direction1[1] ** 2);
                let length2 = Math.sqrt(direction2[0] ** 2 + direction2[1] ** 2);
                // Normalize the direction vector
                direction1 = [direction1[0] / length1, direction1[1] / length1];
                direction2 = [direction2[0] / length2, direction2[1] / length2];
    
                // based on my simple analysis, the direction1 will need to be multiplied in the opposite direction to find the coordinate corresponding to that point. 
                // based on my simple analysis, the direction2 will need to be multiplied in the positive direction to find the coordinate corresponding to that point.
                // Calculate the coordinates of the point at the specified distance from the end
                newPoints.push([surfaceAtMaximumHeight[i][0] - direction1[0] * liftedObjectRadius, surfaceAtMaximumHeight[i][1] - direction1[1] * liftedObjectRadius]);
                newPoints.push([surfaceAtMaximumHeight[i][0] + direction2[0] * liftedObjectRadius, surfaceAtMaximumHeight[i][1] + direction2[1] * liftedObjectRadius]);
                // So at this point, the newPoints should be an array containing the coordinates of the points shifted from the points of the surface, in the order of the points on the surface. 
                // Now you need to take care of the coordinates that need to be identified at the corners because of the circular space at the corners. 
                // Find out the angles of the lines formed by joining the inner point and the point formed by extending by the radial distance. 
                let angles = [];
                angles.push (getAngleFromHorizontal (surfaceAtMaximumHeight[i], newPoints[0]));
                angles.push (getAngleFromHorizontal (surfaceAtMaximumHeight[i], newPoints[1]));
                let angleDifference = 0;
                let adjustedAngle = 0;
                if ((Math.max(...angles)>=1.5*Math.PI && Math.max(...angles)<2*Math.PI) && (Math.min(...angles)>=0 && Math.min(...angles)<0.5*Math.PI)){
                    angleDifference = 2*Math.PI - (Math.max(...angles) - Math.min(...angles));
                    adjustedAngle = Math.max(...angles) + angleDifference/2 ; 
                }
                else {
                    angleDifference = Math.max(...angles) - Math.min(...angles);
                    adjustedAngle = Math.min(...angles) + angleDifference/2 ;
                }
                // let adjustedAngle = Math.min(...angles) + angleDifference/2 ; 
                let x_adjusted = liftedObjectRadius * Math.cos(adjustedAngle);
                let y_adjusted = liftedObjectRadius * Math.sin(adjustedAngle);
                let coordinateOnCircleQuarter = [surfaceAtMaximumHeight[i][0]+x_adjusted, surfaceAtMaximumHeight[i][1]+y_adjusted];
                // Now you need to find the equation of the tangent to the circular sector at that particular point.
                let tangentEquation = getTangentEquation(coordinateOnCircleQuarter, surfaceAtMaximumHeight[i]); 
                let perpendicularEquation1 = getTangentEquation(newPoints[0], surfaceAtMaximumHeight[i]);
                let perpendicularEquation2 = getTangentEquation(newPoints[1], surfaceAtMaximumHeight[i]);
                // Now, you need to solve these equations to find the solutions. 
                newPoints.push(solveLineEquations(tangentEquation, perpendicularEquation1));
                newPoints.push(solveLineEquations(tangentEquation, perpendicularEquation2));
                for (let j=0; j<newPoints.length; j++) {
                    if (j!=1){
                        newCoordinates.push(newPoints[j]);
                    }
                }
                newCoordinates.push(newPoints[1]);
            }
            let edges = [];
            for (let i=0; i<newCoordinates.length-1; i++) {
                //edges.push([newCoordinates[i], newCoordinates[i+1]]);
                edges.push([i, i+1]);
            }
            edges.push([0, newCoordinates.length-1]);
            return {coordinates: newCoordinates,
                    edges: edges};
        }
    
        function findSurfaceWithMaxHeight(polyhedralSurface) {
            let maxSurface = null;
            let maxHeight = -Infinity;
    
            for (let surface of polyhedralSurface) {
                let maxZ = Math.max(...surface.map(vertex => vertex[2]));
                
                if (maxZ > maxHeight) {
                    maxSurface = surface;
                    maxHeight = maxZ;
                }
            }
    
            return maxSurface;
        }
    
        function getCspaceForObstacles(polyhedralSurface, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius) {
            // let adjustedGeometryPolar = getAdjustedGeometryPolar(thetaExtent, adjustedRadius, radius, origin);
            // Finding the top and bottom height of the polyehdral surface. 
            // Initialize minZ and maxZ with the z value of the first vertex
            // So here, it will just return the Cspace object, not the list. 
            let minZ = polyhedralSurface[0][0][2];
            let maxZ = polyhedralSurface[0][0][2];
            // Iterate over each vertex in the polyhedralSurface
            for (var i = 0; i < polyhedralSurface.length; i++) {
            let surface = polyhedralSurface[i];
                for (var j = 0; j < surface.length; j++) {
                    let vertex = surface[j];
                    let z = vertex[2]; // Get the z value of the vertex
                    // Update minZ and maxZ if necessary
                    if (z < minZ) {
                        minZ = z;
                    }
                    if (z > maxZ) {
                        maxZ = z;
                    }
                }
            }
            // Finding the height range of the lifted object. 
            let liftedObjectRange = [];
            liftedObjectRange[0] = setElevation + heigthClearance;
            liftedObjectRange[1] = setElevation + heigthClearance + liftedObjectHeigth;
            // Checking if the object is colliding with the obstacle. 
            let surfaceAtMaximumHeight = [];
            
            if ((liftedObjectRange[0] <= minZ && liftedObjectRange[1] >= minZ) || 
            (liftedObjectRange[0] <= maxZ && liftedObjectRange[1] >= maxZ) ||
            (liftedObjectRange[0] >= minZ && liftedObjectRange[1] <= maxZ)) {
                // Now, you need to get the C-Obstacle for this. 
                // For this, you need the plan coordinates first. 
                surfaceAtMaximumHeight = findSurfaceWithMaxHeight(polyhedralSurface);
                let CSpace = getCSpace(surfaceAtMaximumHeight, liftedObjectRadius);
                return CSpace;
            }
        }
    
        //let setElevation = 6;
        //let heigthClearance = 0.1;
        //let liftedObjectHeigth = 2;
        //let liftedObjectRadius = 1;// It is the radius of rotation of the object itself around its axis. 
        //let Cspaceforobstacles = getCspaceForObstacles(polyhedralSurface, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius);
        //console.log(Cspaceforobstacles)
    
        function ifSurfaceLiesWithinCylinder (radiusCounterweights, surface, origin) {
            //Finding the surface at the top first
            let maxZ = -Infinity;
            let surfaceAtTop = null;
    
            // Iterate over each surface in the polyhedralSurface array
            surface.forEach(surface1 => {
                // Iterate over each vertex in the surface to find the maximum z-coordinate
                surface1.forEach(vertex => {
                    if (vertex[2] > maxZ) {
                        maxZ = vertex[2];
                        surfaceAtTop = surface1;
                    }
                });
            });
            // Check if the distance of each point is less than the radius
            let surfaceLiestWithinCylinder = true; 
            for (let i=0; i<surfaceAtTop.length; i++) {
                let dist = Math.sqrt((surfaceAtTop[i][0]-origin[0])**2 + (surfaceAtTop[i][1]-origin[1])**2);
                if (dist>radiusCounterweights) {
                    surfaceLiestWithinCylinder = false;
                    break;
                }
            }
            return surfaceLiestWithinCylinder;
        }
    
        // Need to check if the counterwights intersect with the existing bounding boxes. The area where intersection is happening needs to be removed. 
        // In this case only the extreme coordinates are required.
        function checkIntersectionOfCounterweights (radiusCounterweights, heightCounterWeights, polyhedralSurfaces, origin, superliftBuffer) {
            // First we get the triangulated representation of the cylinder
            // this function is to return the extent of thetas
            let triangulatedCylinder = distMeshCylinder (origin, radiusCounterweights, heightCounterWeights);
            // console.log(triangulatedCylinder);
            let thetaExtent = [];
            // let intersectionExists = true;
            for (let i=0; i<polyhedralSurfaces.length; i++) {
                let thetaExtentLocal = [];
                // Then we need to find the intersection coordinates of the triangulated cylinder and the polyhedralSurface
                let polyhedralSurface = polyhedralSurfaces[i];
                let coordinatesOfIntersection = detectIntersection (triangulatedCylinder, polyhedralSurface);
                if (coordinatesOfIntersection.intersectionCoordinates.length!=0){
                    // Getting the extreme and innermost coordinates
                    let extremeCoordinates = getExtremeCoordinates(coordinatesOfIntersection.intersectionCoordinates, coordinatesOfIntersection.intersectionCoordinatesEdges, triangulatedCylinder, polyhedralSurface, heightCounterWeights, origin, radiusCounterweights);
                    let innermostCoordinates = getInnermostCoordinates(extremeCoordinates, origin, heightCounterWeights);
                    // Getting the theta extents
                    thetaExtentLocal = getThetaExtents (extremeCoordinates, innermostCoordinates, radiusCounterweights, origin, superliftBuffer);
                    // In the case of the collision of the counterweights, the theta extents need to be adjusted to the opposite side because the intersection of the counterweight will resist the boom operation. 
                    // So, in every theta extent, I add pie. 
                    for (let i=0; i<thetaExtentLocal.length; i++) {
                        thetaExtentLocal[i] = thetaExtentLocal[i] + Math.PI;
                        if (thetaExtentLocal[i]>2*Math.PI) {
                            thetaExtentLocal[i] = thetaExtentLocal[i] - 2*Math.PI;
                        }
                    }
                    thetaExtent.push(thetaExtentLocal);
                }
                else if (ifSurfaceLiesWithinCylinder(radiusCounterweights, surface, origin)=true){
                    // This is the case when the whole bounding box actually lies within the cylindrical volume
                    // So, in this case, the theta extents have to directly identified from the object itself. 
                    // So the theta extents will be basically the points with minimum and maximum angle of rotation from the centre. 
                    // If the nearest points lie in the fourth and first quadrant, it will be slightly different. 
                    thetaExtentLocal = getThetaExtentsCounterweight (surface, radiusCounterweights, origin, superliftBuffer);
                    for (let i=0; i<thetaExtentLocal.length; i++) {
                        thetaExtentLocal[i] = thetaExtentLocal[i] + Math.PI;
                        if (thetaExtentLocal[i]>2*Math.PI) {
                            thetaExtentLocal[i] = thetaExtentLocal[i] - 2*Math.PI;
                        }
                    }
                    thetaExtent.push(thetaExtentLocal);
                }
                /*else {
                    intersectionExists = false;
                    thetaExtent.push(0);
                }*/
            }
            return thetaExtent;
        }
    
        //var heightCrane = 5;
        //let radiusCounterweights = 4;
        //let heightCounterWeights = 3;
        //let superliftBuffer = 0.5;
        //let coordinatesOfIntersection = detectIntersection (triangulatedCylinder, polyhedralSurfaces[0]);
        // console.log(coordinatesOfIntersection.intersectionCoordinates);
        //let extremeCoordinates = getExtremeCoordinates(coordinatesOfIntersection.intersectionCoordinates, coordinatesOfIntersection.intersectionCoordinatesEdges, triangulatedCylinder, polyhedralSurfaces[0], heightCounterWeights, origin);
        //console.log(extremeCoordinates.extremeCoordinates);
        //let innermostCoordinates = getInnermostCoordinates(extremeCoordinates, origin, heightCounterWeights);
        //console.log(innermostCoordinates);
    
        //let counterweightsIntersection = checkIntersectionOfCounterweights (radiusCounterweights, heightCounterWeights, polyhedralSurfaces, origin, superliftBuffer);
        //console.log(counterweightsIntersection);
        // Finding out the updated geometry for constraint checking
    
        // What you need to do is to run all the above code for every polyhedral surface in the model. 
        // The Theta extent and adjusted radius should be for all the bounding boxes. 
        // In this way, you can divide the whole circle into a set of sectors. And create the edges for each sector individually. 
        // For this, the thetaExtent and adjustedRadius should contain all the values, instead of only one value. You need to run them on some sort of loop. 
        // One important thing here is that the thetaExtent just stores the extents of theta for which the radius is adjusted radius. 
        /*let polyhedralSurfaces =  [[[[-2, 2, 4],[-2, 6, 4],[2, 6, 4], [2, 2, 4], [-2, 2, 4]],
                                    [[-2, 2, 4], [-2, 2, 0], [2, 2, 0], [2, 2, 4], [-2, 2, 4]],
                                    [[-2, 6, 4], [-2, 6, 0], [2, 6, 0], [2, 6, 4], [-2, 6, 4]],
                                    [[-2, 2, 4], [-2, 2, 0], [-2, 6, 0], [-2, 6, 4], [-2, 2, 4]],
                                    [[2, 2, 4], [2, 6, 4], [2, 6, 0], [2, 2, 0], [2, 2, 4]]] , 
                                [[[1, -1, 8],[1, -6, 8],[6, -6, 8], [6, -1, 8], [1, -1, 8]],
                                    [[1, -1, 8], [1, -1, 0], [6, -1, 0], [6, -1, 8], [1, -1, 8]],
                                    [[1, -6, 8], [1, -6, 0], [6, -6, 0], [6, -6, 8], [1, -6, 8]],
                                    [[1, -1, 8], [1, -1, 0], [1, -6, 0], [1, -6, 8], [1, -1, 8]],
                                    [[6, -1, 8], [6, -6, 8], [6, -6, 0], [6, -1, 0], [6, -1, 8]]]];
        let origin = [0,0,0];
        let radiusOfCone = 5;
        let heightOfCone = 5;
        let boomClearance = 0.5;
        let boomLength = 7.07;
        let setElevation = 3;
        let heigthClearance = 0.1;
        let liftedObjectHeigth = 1;
        let liftedObjectRadius = 2;
        let radiusCounterweights = 3;
        let heightCounterWeights = 1.5;
        let superliftBuffer = 0.5;*/
    
        function getSpaceAdjustments (polyhedralSurfaces, origin, radiusOfCone, heightOfCone, boomClearance, boomLength, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius,
                                            radiusCounterweights, heightCounterWeights, superliftBuffer) {
            // So you have to iterate over all the obstructions to perform the checking.
            // What you can do is, instead of checking for every obstruction, you can put a simple check if the obstruction lies within the range of movement of the crane.
            let triangulatedSurface = distMeshCone (origin, radiusOfCone, heightOfCone);
            /*for (let i=0; i<triangulatedSurface.finalPoints.length; i++) {
                console.log(triangulatedSurface.finalPoints[i]);
            } 
            for (let i=0; i<triangulatedSurface.bars.length; i++) {
                console.log(triangulatedSurface.bars[i]);
            }*/
            let thetaExtent = [];
            let adjustedRadius = [];
            let Cspace = [];
            for (let i=0; i<polyhedralSurfaces.length; i++) {
                let surface = polyhedralSurfaces[i];
                // Getting the intersection coordinates
                let coordinatesOfIntersection = detectIntersection (triangulatedSurface, surface);
                if (coordinatesOfIntersection.intersectionCoordinates.length!=0) {
                    // Getting the extremeCoordinates
                    let extremeCoordinates = getExtremeCoordinates(coordinatesOfIntersection.intersectionCoordinates, coordinatesOfIntersection.intersectionCoordinatesEdges, triangulatedSurface, surface, heightOfCone, origin, radiusOfCone);
                    // Getting the innermostCoordinates
                    let innermostCoordinates = getInnermostCoordinates(extremeCoordinates, origin, heightOfCone);
                    thetaExtent.push(getThetaExtents (extremeCoordinates, innermostCoordinates, radiusOfCone, origin, boomClearance));
                    let heightOfPolyhedralSurface = findSurfaceWithMaxHeight(surface)[0][2];
                    adjustedRadius.push(getAdjustedRadius (boomLength, innermostCoordinates, origin, boomClearance, heightOfPolyhedralSurface, heightOfCone));
                }
                Cspace.push(getCspaceForObstacles (surface, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius));
            }
            let radiusAngleRanges = getAdjustedGeometryPolar (thetaExtent, adjustedRadius, radiusOfCone, origin);
            let thetaExtentCounterweights = checkIntersectionOfCounterweights (radiusCounterweights, heightCounterWeights, polyhedralSurfaces, origin, superliftBuffer);
            
            return {radiusAngleRanges : radiusAngleRanges,
                    Cspace : Cspace, 
                    thetaExtentCounterweights : thetaExtentCounterweights};
            // This way, this function should return the overall circular geometry, with different sectors and their radius values and the Cspace. 
            // radiusAngleRanges is a list, where each element is a list. The contents of each list are [theta1, theta2, radius]. 
            // Cspace will be a list of objects, where each object contains a list of coordinates of points and the pairs of coordinates connecting edges. 
            // Wherever you are using it, you need to check if the CSpace actually contains some value or not. 
            // This is because the Cspace is relevant only of the obstacle height is greater than the set elevation. 
        }
    
        let adjustedSpace = getSpaceAdjustments (polyhedralSurfaces, origin, radiusOfCone, heightOfCone, boomClearance, boomLength, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius,
            radiusCounterweights, heightCounterWeights, superliftBuffer);
        // console.log(adjustedSpace);
    
        function checkIfPointLiesWithinRadiusAngleRange (pointAngle, radiusAngleRange, point, origin) {
            // Here, you need to first find out in which sector the point lies and then check if the point is at a distance less than the corresponding radius of that sector. 
            // If it is, then return true, otherwise, return false. 
            for (let i=0; i<radiusAngleRange.length; i++) {
                let angles = [radiusAngleRange[i][0], radiusAngleRange[i][1]];
                if ((Math.max(...angles)>=1.5*Math.PI && Math.max(...angles)<2*Math.PI) &&
                (Math.min(...angles)>=0 && Math.max(...angles)<0.5*Math.PI)) {
                    if ((pointAngle>=0 && pointAngle<Math.min(...angles)) || 
                    (pointAngle>=Math.max(...angles) && pointAngle<2*Math.PI)){
                        if (Math.sqrt((point[0]-origin[0])**2 + (point[1]-origin[1])**2) <= radiusAngleRange[i][2]) {
                            return true;
                        }
                        else {
                            return false;
                        }
                    }
                }
                else {
                    if (pointAngle>=Math.min(...angles) && pointAngle<Math.max(...angles)){
                        if (Math.sqrt((point[0]-origin[0])**2 + (point[1]-origin[1])**2) <= radiusAngleRange[i][2]) {
                            return true;
                        }
                        else {
                            return false;
                        }
                    }
                }
            }
        }
    
        function checkIfPointLiesWithinAngleRange (pointAngle, thetaExtentsRemoved) {
            // here, you are just trying to check if the point lies within a range of thetas, which are removed from the circular geometry. 
            for (let i=0; i<thetaExtentsRemoved.length; i++) {
                let angles = [thetaExtentsRemoved[i][0], thetaExtentsRemoved[i][1]];
                if ((Math.max(...angles)>=1.5*Math.PI && Math.max(...angles)<2*Math.PI) &&
                (Math.min(...angles)>=0 && Math.min(...angles)<0.5*Math.PI)) {
                    if ((pointAngle>=0 && pointAngle<Math.min(...angles)) || 
                    (pointAngle>=Math.max(...angles) && pointAngle<2*Math.PI)){
                            return true;
                        }
                        // else {
                        //    return false;
                        //}
                    }
                else {
                    if (pointAngle>=Math.min(...angles) && pointAngle<Math.max(...angles)){
                            return true;
                        }
                        //else {
                        //    return false;
                    //}
                }
            }
            return false;
        }
    
        function checkIfPointLiesWithinCSpace (point, Cspace) {
            // So here, the objective is to check if the point lies within the Cspace. 
            // CSpace has a set of geometries. Each object in the CSpace list should contain a list of coordinates and a list of lists containing edges.
            // A simple thing you can check here is, if the ray originating from the point intersects odd or even number of edges in the CSpace polygons.  
            for (let i=0; i<Cspace.length; i++) {
                let ifPointLiesWithinCspace = checkIfLiesWithinLimit (Cspace[i], point, "polygon");
                if (ifPointLiesWithinCspace == true) {
                    return true;
                }
            }
            return false;
        }
    
        function checkIfRotationPossible (pickPoint, setPoint, thetaExtentsRemoved, origin) {
            // This function is to check if the rotation is possible based on the pickpoints, set points, removed theta extents and the origin of rotation. 
            // Getting the pick and set angles
            let pickPointAngle = getAngleFromHorizontal(origin, pickPoint);
            let setPointAngle = getAngleFromHorizontal(origin, setPoint);
            // There will always be two angle ranges betweent the two thetas. 
            // One, which will go from the minimum to maximum. 
            // Second, which will go from the Maximum to zero and then zero to minimum. 
            
            let emptyRangeLiesWithinPickSetPoints = 0;
    
            // Checking for the first range. 
            for (let i=0; i<thetaExtentsRemoved.length; i++) {
                if ((Math.min(pickPointAngle, setPointAngle) <= Math.min(...thetaExtentsRemoved[i])) &&
                (Math.max(pickPointAngle, setPointAngle) >= Math.max(...thetaExtentsRemoved[i]))) {
                    emptyRangeLiesWithinPickSetPoints+=1;
                    break;
                    // As soon as you get a sector in between, you break from the loop. 
                }
            }
    
            // Checking for the second range. 
            // Here you will have to go from the maximum theta to zero and then zero to minimum theta. 
            let range1 = [Math.max(pickPointAngle, setPointAngle), 2*Math.PI];
            let range2 = [0, Math.min(pickPointAngle, setPointAngle)];
            for (let i=0; i<thetaExtentsRemoved.length; i++)  {
                // Three cases will be here.
                // First, both angles lie in range 1
                if ((range1[0] <= Math.min(...thetaExtentsRemoved[i])) &&
                (range1[1] >= Math.max(...thetaExtentsRemoved[i]))) {
                    emptyRangeLiesWithinPickSetPoints+=1;
                    break;
                }
                // Second, both angles lie in range 2
            
                else if ((range2[0] <= Math.min(...thetaExtentsRemoved[i])) &&
                (range2[1] >= Math.max(...thetaExtentsRemoved[i]))) {
                    emptyRangeLiesWithinPickSetPoints+=1;
                    break;
                }
    
                // Third, one angle lie in range 1 and second in range2. 
                else if ((Math.min(pickPointAngle, setPointAngle) >= Math.min(...thetaExtentsRemoved[i])) && 
                (Math.min(...thetaExtentsRemoved[i])>=0) && 
                (Math.max(pickPointAngle, setPointAngle) <= Math.max(...thetaExtentsRemoved[i])) &&
                (Math.max(...thetaExtentsRemoved[i]) <=2*Math.PI)) {
                    emptyRangeLiesWithinPickSetPoints+=1;
                    break;
                }
            }
    
            if (emptyRangeLiesWithinPickSetPoints == 2) {
                return false;
            }
            else {
                return true;
            }
        }
    
        /*function getAdjustedGeometry (thetaExtent, adjustedRadius, radius, origin) {
            // This function should give you the edges for all the sectors in the circle.
    
            // Sorting the theta extents based on the first value of of each inner list.
            let extentTheta =  thetaExtent;
            let sortedThetaExtent = extentTheta.sort((a, b) => a[0] - b[0]);
            let radiusAngleRanges = [];
            // So here, for all the values in the sortedThetaExtent, you store the radius along with it. 
            // First storing for those values which have the radius less than the outer radius.
            for (let i=0; i<sortedThetaExtent.length; i++) {
                for (let j=0; j<thetaExtent.length; j++) {
                    if (JSON.stringify(sortedThetaExtent[i]) === JSON.stringify(thetaExtent[j])) {
                        radiusAngleRanges.push([sortedThetaExtent[i][0], sortedThetaExtent[i][1], adjustedRadius[j]]);
                    }
                }
            }
            // Then storing the values for which the radius is equal to the outer radius. 
            for (let i=0; i<sortedThetaExtent.length-1; i++) {
                radiusAngleRanges.push([sortedThetaExtent[i][1], sortedThetaExtent[i+1][0], radius]);
            }
            // Adding the value corresponding to the last sortedThetaExtent and the first one. 
            radiusAngleRanges.push([sortedThetaExtent[sortedThetaExtent.length-1][1], sortedThetaExtent[0][0], radius]);
    
            // Now, you need to create the lines for every sector. 
            let edgeRepresentationOfArcs = [];
            for (let i=0; i<radiusAngleRanges.length; i++) {
                edgeRepresentationOfArcs.push(getedgeRepresentationOfArcs(radiusAngleRanges[i]));
            }
            let additionalEdges = [];
            for (let i=0; i<edgeRepresentationOfArcs.length; i++) {
                additionalEdges.push([edgeRepresentationOfArcs[i][edgeRepresentationOfArcs[i].length-1][1], edgeRepresentationOfArcs[i+1][0][0]]);
            }
            edgeRepresentationOfArcs.push(additionalEdges)
            return edgeRepresentationOfArcs;
        }*/
    
        function checkIfLiftIsFeasible (polyhedralSurfaces, origin, radiusOfCone, heightOfCone, boomClearance, boomLength, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius,
            radiusCounterweights, heightCounterWeights, superliftBuffer, pickPoint, setPoint, minimumRadiusCone) {
            // This function is to finally check if the path can actually exist between the pickArea and setArea or not. 
            // Getting the feasible operating range first. 
            let feasibleOperatingSpace = getSpaceAdjustments (polyhedralSurfaces, origin, radiusOfCone, heightOfCone, boomClearance, boomLength, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius,
                radiusCounterweights, heightCounterWeights, superliftBuffer);
            let radiusAngleRange = feasibleOperatingSpace.radiusAngleRanges;
            let Cspace = feasibleOperatingSpace.Cspace;
            let thetaExtentCounterweights = feasibleOperatingSpace.thetaExtentCounterweights;
            let thetaExtentsRemoved = thetaExtentCounterweights;
            // In the radiusAngleRange, if any radius is less than Rmin, the operation is not possible within that range. 
            // One thing is missing here. That is of the C-Space go beyond 
            for (let i=0; i<radiusAngleRange; i++) {
                if (radiusAngleRange[i][2]<=minimumRadiusCone) {
                    thetaExtentsRemoved.push(radiusAngleRange[i][0], radiusAngleRange[i][1]);
                }
            }
    
            
            let pickPointAngle = getAngleFromHorizontal(origin, pickPoint);
            let setPointAngle = getAngleFromHorizontal(origin, setPoint);
            
            // We have to check if these points lie within the removed sectors or not. 
            let pickPointWithinRemovedSectorLimits = checkIfPointLiesWithinAngleRange (pickPointAngle, thetaExtentsRemoved); // The value will be true if the point lies within the range.
            let setPointWithinRemovedSectorLimits = checkIfPointLiesWithinAngleRange (setPointAngle, thetaExtentsRemoved); // The value will be true if the point lies within the range.
            
    
            // If the point lies within the removed sectors, we do not need to check any further. We can just return from here. 
            if (pickPointWithinRemovedSectorLimits == true || setPointWithinRemovedSectorLimits== true) {
                return false;
            }
    
            // Next, you need to check if the pick point and set point lie within the C-Spaces. 
            // Remember that there are multiple C-Spaces in the Cspace. 
            let pickPointWithinCspace = checkIfPointLiesWithinCSpace (pickPoint, Cspace); // The value of this variable will be true if the point lies within the Cspace. 
            let setPointWithinCspace = checkIfPointLiesWithinCSpace (setPoint, Cspace); // The value of this variable will be true if the point lies within the Cspace.
            // If the point lies within the Cspaces, we do not need to check any further. We can just return from here. 
            if (pickPointWithinCspace == true || setPointWithinCspace== true) {
                return false;
            }
            
            // If the above two conditions are satisfied, we check if the pick and set points lie within the ranges of the radiusAngles.
            let pickPointWithinLimits = checkIfPointLiesWithinRadiusAngleRange (pickPointAngle, radiusAngleRange, pickPoint, origin); // The value will be true if the point lies within the range. 
            let setPointWithinLimits = checkIfPointLiesWithinRadiusAngleRange (setPointAngle, radiusAngleRange, setPoint, origin); // The value will be true if the point lies within the range. 
            // pickPointWithinLimits and setPointWithinLimits stores if the point lies within the radiusAngleRange. 
            // If the point lies within the Cspaces, we do not need to check any further. We can just return from here. 
            if (pickPointWithinLimits == false || setPointWithinLimits == false) {
                return false;
            }
    
            // Finally, if all the above conditions are met, we go further to check the path can exist between the two points or not. 
            // For the path, we actually do not have anything here. 
            // We need to check if between the pickpoint and set point, there exists empty areas such that there is no possibility for the crane
            // to perform the operation. 
            // What we need to check is that if the empty areas exist on the both sides between the pick and set points. 
            // This means that there is no space between the two, between which we can rotate the crane. 
            let rotationPossible = checkIfRotationPossible (pickPoint, setPoint, thetaExtentsRemoved, origin);
            if (rotationPossible == false) {
                return false;
            }
            else if (rotationPossible == true) {
                return true;
            }
        }
        // This function is to check if there exisits a path between the pick point and the set point. 
    
        //let setPoint = [-2, 1];
        //console.log(getAngleFromHorizontal(origin, setPoint));
        //let pickPoint = [1.5, -2];
        //console.log(getAngleFromHorizontal(origin, pickPoint));
        // let origin = [0,0,0];
        //let thetaExtentsRemoved = [[4*Math.PI/3, 5*Math.PI/3], [7*Math.PI/4, Math.PI/8]];
        //let ifRotationPossible = checkIfRotationPossible(pickPoint, setPoint, thetaExtentsRemoved, origin);
        //console.log(ifRotationPossible);
        // let point1 = [2, 2, 0];
        // let point2 = [-3, 5, 0];
        // let closestPoint = getClosestPoint(point1, point2, origin);
        // console.log(closestPoint);
        //let thetaextents = getThetaExtents (extremeCoordinates, innermostCoordinates, radius, origin, boomClearance);
        //console.log(thetaextents);
        return (checkIfLiftIsFeasible (polyhedralSurfaces, origin, radiusOfCone, heightOfCone, boomClearance, boomLength, setElevation, heigthClearance, liftedObjectHeigth, liftedObjectRadius,
            radiusCounterweights, heightCounterWeights, superliftBuffer, pickPoint, setPoint, minimumRadiusCone));
    })
    .catch(error => {
        console.error('Failed to load the d3 module:', error);
    }); 
}

// Function for detecting intersection between the oriented bounding boxes. 
function detectIntersectionOBB (obb1, obb2) {
    class Vector3 {
        constructor(x, y, z) {
            this.x = x;
            this.y = y;
            this.z = z;
        }
    
        subtract(vector) {
            return new Vector3(this.x - vector.x, this.y - vector.y, this.z - vector.z);
        }
    
        dot(vector) {
            return this.x * vector.x + this.y * vector.y + this.z * vector.z;
        }
    
        cross(vector) {
            return new Vector3(
                this.y * vector.z - this.z * vector.y,
                this.z * vector.x - this.x * vector.z,
                this.x * vector.y - this.y * vector.x
            );
        }
    
        normalize() {
            const length = Math.sqrt(this.x * this.x + this.y * this.y + this.z * this.z);
            return new Vector3(this.x / length, this.y / length, this.z / length);
        }
    }
    
    class OBB {
        constructor(center, halfWidths, axes) {
            this.center = center; // Vector3
            this.halfWidths = halfWidths; // Vector3
            this.axes = axes; // Array of Vector3, representing the local x, y, z axes
        }
    
        getProjection(axis) { // This is finding the projection of the given box along a particular axis, which is an argument to this function.
            let radius = this.halfWidths.x * Math.abs(axis.dot(this.axes[0])) +
                         this.halfWidths.y * Math.abs(axis.dot(this.axes[1])) +
                         this.halfWidths.z * Math.abs(axis.dot(this.axes[2]));
            let projectionCenter = this.center.dot(axis); // This will find the projection of the centre of the box.
            return radius
            //{ min: projectionCenter - radius, max: projectionCenter + radius };
        }
    }
    
    function getProjection2(coordinate, axis) {
        return coordinate.dot(axis)
    }
    
    function doOBBSIntersect(obbA, obbB) {
        let axes = [
            ...obbA.axes,
            ...obbB.axes,
            ...obbA.axes.map(axisA => obbB.axes.map(axisB => axisA.cross(axisB))).flat()
        ];
    
        for (let axis of axes) {
            if (axis.x === 0 && axis.y === 0 && axis.z === 0) continue; // Skip zero length axis from parallel cross products
            axis = axis.normalize(); // Normalize the axis for accurate projections
    
            let projectionA = obbA.getProjection(axis);
            let projectionB = obbB.getProjection(axis);
            let projectionCenter = getProjection2(obbB.center.subtract(obbA.center), axis)
    
            if (projectionCenter > projectionA + projectionB) {
                return false;
            }
            //if (projectionA.max < projectionB.min || projectionB.max < projectionA.min) {
            //   return false; // No intersection found on this axis
            //}
        }
        return true; // All axes tested, no separating axis found
    }
    
    // Example usage
    const obbA = new OBB(new Vector3(0, 0, 0), new Vector3(1, 1, 1), [
        new Vector3(1, 0, 0),
        new Vector3(0, 1, 0),
        new Vector3(0, 0, 1)
    ]); // So this creates an object obbA of class OBB, which takes the arguments as centre - new Vector3(0, 0, 0), halfwidths - new Vector3(1, 1, 1), and axes - an array of three Vector3 objects
    // Thet define the local coordinate system of the bounding box. These vectors represent the directions of bounding box's local X, U, Z axes. 
    
    const obbB = new OBB(new Vector3(0, 0, 0), new Vector3(1, 1, 1), [
        new Vector3(1, 0, 0),
        new Vector3(0, 1, 0),
        new Vector3(0, 0, 1)
    ]);
    
    var intersection = doOBBSIntersect(obbA, obbB);
    console.log(intersection);
    return intersection;
    //console.log(doOBBSIntersect(obbA, obbB)); // Output: true or false based on intersection
    //return 
}

