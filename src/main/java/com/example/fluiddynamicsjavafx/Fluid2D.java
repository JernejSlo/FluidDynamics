package com.example.fluiddynamicsjavafx;

import javafx.animation.AnimationTimer;
import javafx.application.Application;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.scene.shape.Line;
import javafx.stage.Stage;

import java.io.IOException;
import java.util.*;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;
import java.util.concurrent.Future;
import java.util.concurrent.TimeUnit;


public class Fluid2D extends Application {

    // Project by : Jernej Koprivnikar
    //
    // Credit to the creator of this video - https://www.youtube.com/watch?v=-0m05gzk8nk&ab_channel=MachineLearning%26Simulation
    // I translated the code into java and added surface tension and foam. Original code is in python and mostly uses python
    // libraries to do the calculations, so this code is not that similar. The things that are the same are the equations which are
    // from this paper: https://matthias-research.github.io/pages/publications/sca03.pdf


    // CHANGEABLE PARAMETERS

    public static int SCREEN_WIDTH = 800;
    public static int SCREEN_HEIGHT = 600;

    public boolean paralel = false;

    public boolean emitter = true;
    public int[] emitterPosition = {SCREEN_WIDTH/2,SCREEN_HEIGHT/2};
    public int waitBetweenEmits = 10;
    public int emmitParticlesNum = 10;

    public static boolean turnOnSurfaceTension = false;
    public static boolean colorfull = true;
    public static boolean decrease_foam_volume = true;
    static int dam = 0;
    static boolean simulateDamBreak = false;

    public static int PARTICLE_RADIUS = 8;
    public static int NUM_PARTICLES = (SCREEN_WIDTH+SCREEN_HEIGHT)/PARTICLE_RADIUS*5;
    public static int num_threads = 30;
    protected static ArrayList<KdTree.Node> [] neighbors = new ArrayList[NUM_PARTICLES];
    double[] startVector = {500,-200};


    // SOME CAN BE CHANGED BUT REALLY SHOULDN'T BE
    public static ExecutorService executorService;

    public static KdTree tree;

    private static List<ParticleDrawn> particles = new ArrayList<>();
    private List<Line> curves = new ArrayList<>();
    public static double[][] newForces = new double[NUM_PARTICLES][2];

    public static double [][] positions = new double[NUM_PARTICLES][2];
    public static List<KdTree.Node> particleCoordinates = new ArrayList<>();

    static double smoothingLength = PARTICLE_RADIUS*2;
    Pane root;


    static double mass = 1;
    double [] zeros;
    static double timeStep = 0.01;

    public static double ISOTROPIC_EXPONENT = 20;
    public static double BASE_DENSITY = 1;
    public static double [] CONSTANT_FORCE = {0,-0.1};
    public static double DYNAMIC_VISCOSITY = 0.5;

    public static double NORMALIZATION_DENSITY = (315 * mass) / (64 * Math.PI * Math.pow((smoothingLength),9));
    public static double NORMALIZATION_VISCOUS_FORCE = (45 * DYNAMIC_VISCOSITY * mass) / (Math.PI * (Math.pow((smoothingLength),6)));
    public static double NORMALIZATION_PRESSURE_FORCE = -((45 * mass) / (Math.PI * Math.pow((smoothingLength),6)));

    public static final int MIN_X = 0;
    private static final int MIN_Y = 0;
    public static double MAX_X = SCREEN_WIDTH;
    protected static double MAX_Y = SCREEN_HEIGHT;

    // The coefficient of restitution for the vectors
    public static final double COR = 0.9;

    @Override
    public void start(Stage stage) throws IOException {
        Scene scene = new Scene(createContent());

        stage.setTitle("Fluid Dynamics");



        stage.setScene(scene);
        stage.show();


    }

    private Parent createContent(){

        zeros = new double[2];
        zeros[0] = 0;
        zeros[1] = 0;

        Pane root = new Pane();
        root.setPrefSize(SCREEN_WIDTH,SCREEN_HEIGHT);
        root.setStyle("-fx-background-color: black;");

        for (int i = 0; i < num_threads; i++) {
            var curve = new Line();
            curves.add(curve);
            root.getChildren().add(curve);
        }
        this.root = root;
        if (!emitter){
            for (int i = 0; i < NUM_PARTICLES; i++) {
                var p = new ParticleDrawn();
                particles.add(p);
                root.getChildren().add(p);
            }


            double x = 0;
            double y = 0;
            if (simulateDamBreak){
                for (int i = 0; i < NUM_PARTICLES; i++) {
                    ParticleDrawn particle = particles.get(i);
                    double num = smoothingLength;

                    particle.setTranslateX(x);
                    particle.setTranslateY(y);
                    positions[i][0] = x; //+Math.random()*1;
                    positions[i][1] = y;

                    if (x >= (MAX_X)/2){
                        x = 0;
                        y += num;
                    }

                    x+= num/2;

                }}
            else {
                x = MAX_X/4;
                y = MAX_Y/4;
                for (int i = 0; i < NUM_PARTICLES; i++) {
                    ParticleDrawn particle = particles.get(i);
                    double num = smoothingLength;

                    particle.setTranslateX(x);
                    particle.setTranslateY(y);
                    positions[i][0] = x; //+Math.random()*1;
                    positions[i][1] = y;

                    if (x >= (3*MAX_X)/4){
                        x = MAX_X/4;
                        y += num;
                    }

                    x+= num/2;

                }
            }
            for (int i = 0; i < positions.length; i++) {
                KdTree.Node newNode = new KdTree.Node(new double[] {positions[i][0],positions[i][1]});
                double [] velocity = startVector;

                newNode.setVelocity(velocity);
                newNode.setDensity(BASE_DENSITY);
                newNode.id = i;
                particleCoordinates.add(newNode);

            }
        }
        else {
            Random rand = new Random();
            for (int i = 0; i < emmitParticlesNum; i++) {
                var p = new ParticleDrawn();
                particles.add(p);
                root.getChildren().add(p);
            }

            int x = emitterPosition[0];
            int y = emitterPosition[1];
            for (int i = 0; i < emmitParticlesNum; i++) {
                ParticleDrawn particle = particles.get(i);
                double num = smoothingLength;

                particle.setTranslateX(x);
                particle.setTranslateY(y);
                positions[i][0] = x; //+Math.random()*1;
                positions[i][1] = y;



                x+= num/2;

            }

            for (int i = 0; i < emmitParticlesNum; i++) {
                KdTree.Node newNode = new KdTree.Node(new double[] {positions[i][0],positions[i][1]});
                double [] velocity = {rand.nextInt(200)+1,rand.nextInt(200)+1};

                newNode.setVelocity(velocity);
                newNode.setDensity(BASE_DENSITY);
                newNode.id = i;
                particleCoordinates.add(newNode);

            }
        }




        tree = new KdTree(2,particleCoordinates);
        tree.setRadius(smoothingLength);



        AnimationTimer timer = new AnimationTimer() {
            @Override
            public void handle(long now) {
                if (paralel){
                    onUpdateParalel();
                }
                else {

                    onUpdate();
                }

            }
        };
        timer.start();

        return root;
    }

    private static class ParticleDrawn extends Parent {
        Circle shape;

        ParticleDrawn(){
            shape = new Circle(PARTICLE_RADIUS*2);
            shape.setFill(Color.LIGHTSKYBLUE);

            getChildren().addAll(shape);
        }
        public void setColor(Color color){
            shape.setFill(color);
        }

        public void setRadius(int radius){
            shape.setRadius(radius);
        }
    }



    public static void calculateDensity(KdTree.Node currParticle){
        ArrayList collisions = neighbors[currParticle.id];
        double [] pl = currParticle.getCoords_();
        double density = NORMALIZATION_DENSITY * Math.pow((Math.pow(smoothingLength,2)-Math.pow(0,2)),3);
        for (int i = 0; i < collisions.size(); i++) {
            KdTree.Node neighbour = (KdTree.Node) collisions.get(i);
            double [] c = neighbour.getCoords_();
            double distance = getDistance(pl[0],pl[1],c[0],c[1]);

            density+=NORMALIZATION_DENSITY * Math.pow((Math.pow(smoothingLength,2)-Math.pow(distance,2)),3);

        }
        currParticle.setDensity(density);
    }

    public static void calculatePressure(KdTree.Node currParticle){
        double pressure = ISOTROPIC_EXPONENT * (currParticle.getDensity() -BASE_DENSITY);
        currParticle.setPressure(pressure);
    }

    public double[] normalize(double[] vector){
        double magnitude = Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1]);
        double[] normalizedVector = {vector[0]/magnitude, vector[1]/magnitude};



        return normalizedVector;
    }
    public static double [] sub(double[] vector, double[] subtract){
        double [] subbedVector = new double[2];

        subbedVector[0] = vector[0]-subtract[0];
        subbedVector[1] = vector[1]-subtract[1];

        return subbedVector;
    }
    public static double [] add(double[] vector, double[] add){
        double [] addedVector = new double[2];

        addedVector[0] = vector[0]+add[0];
        addedVector[1] = vector[1]+add[1];

        return addedVector;
    }
    public static double [] mul(double[] vector1, double scalar){
        double [] mulVector = new double[2];

        mulVector[0] = vector1[0]*scalar;
        mulVector[1] = vector1[1]*scalar;

        return mulVector;
    }
    public double [] vecmul(double [] vector1,double [] vector2){
        double [] mulVector = new double[2];

        mulVector[0] = vector1[0]*vector2[0];
        mulVector[1] = vector1[1]*vector2[1];

        return mulVector;
    }
    public double [] pow(double [] vector1,double pow){
        double [] mulVector = new double[2];

        mulVector[0] = Math.pow(vector1[0],pow);
        mulVector[1] = Math.pow(vector1[1],pow);

        return mulVector;
    }
    public static double [] div(double[] vector1, double scalar){
        double [] divVector = new double[2];

        divVector[0] = vector1[0]/scalar;
        divVector[1] = vector1[1]/scalar;

        return divVector;
    }
    public double dot(double [] vector1,double [] vector2){

        double dotProduct = (vector1[0]* vector2[0]) + (vector1[1]* vector2[1]);

        return dotProduct;
    }

    public static void calculateForces(KdTree.Node currParticle) {
        double[] velocity = {0,0};
        double[] pressureForce = {0,0};
        double[] viscousForce = {0,0};
        double[] currentVelocity = currParticle.getVelocity();
        ArrayList collisions = neighbors[currParticle.id];
        double [] pl = currParticle.getCoords_();
        for (int i = 0; i < collisions.size(); i++) {
            KdTree.Node neighbour = (KdTree.Node) collisions.get(i);
            double [] c = neighbour.getCoords_();
            double distance = getDistance(pl[0],pl[1],c[0],c[1]);
            double scalar = (neighbour.getPressure()+currParticle.getPressure() ) / (2* neighbour.getDensity()) * (Math.pow((smoothingLength - distance),2));
            double [] forceVectorPressure = div(sub(c,pl),distance);
            pressureForce = sub(velocity,
                    (mul(
                            mul(forceVectorPressure,scalar)
                            ,NORMALIZATION_PRESSURE_FORCE))
            );

            //adding viscous force
            forceVectorPressure = div(sub(neighbour.getVelocity(),currentVelocity),neighbour.getDensity());
            scalar = smoothingLength-distance;

            viscousForce = add(velocity,
                    (mul(
                            mul(forceVectorPressure,scalar)
                            ,NORMALIZATION_VISCOUS_FORCE))
            );

        }


        velocity = add(add(pressureForce,viscousForce),CONSTANT_FORCE);


        double [] velocityCalculated = add(currentVelocity, div(mul(velocity,timeStep),currParticle.getDensity()));
        newForces[currParticle.id] = velocityCalculated;
        //currParticle.setVelocity(velocityCalculated);
        //currParticle.setCoords_(add(pl,mul(currParticle.getVelocity(),timeStep)));

    }

    public static void damBreak(KdTree.Node currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] >= MAX_X/2+smoothingLength){

            double [] newVelocity = {-velocity[0]*COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {MAX_X/2+smoothingLength, currparticle.getCoords_()[1]});
        }
    }

    // dodaj da se odbije od zida
    public static void checkIfBounced(KdTree.Node currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] <= MIN_X){

            double [] newVelocity = {-velocity[0]*COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {0, currparticle.getCoords_()[1]});
        }
        if (coords[0] >= MAX_X){

            double [] newVelocity = {-velocity[0]*COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {MAX_X, currparticle.getCoords_()[1]});
        }
        if (coords[1] <= 0){
            //double [] newVelocity = {velocity[0],-velocity[1]};
            double [] newVelocity = {velocity[0],-velocity[1]*COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],0});
        }
        if (coords[1] >= MAX_Y){

            double [] newVelocity = {velocity[0],-velocity[1]*COR};
            //double [] newVelocity = {velocity[0],-velocity[1]};
            currparticle.setVelocity(newVelocity);

            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],MAX_Y});
        }
    }

    public void updateParticles() {
        for (KdTree.Node p : particleCoordinates) {



            particles.get(p.id).setTranslateX(p.getCoords_()[0]);
            particles.get(p.id).setTranslateY(Math.abs(p.getCoords_()[1]-MAX_Y));
            particles.get(p.id).setColor(p.getColor());
            particles.get(p.id).setRadius(p.getRadius());
            //particles.get(p.id).setColor(Color.TRANSPARENT);
            //particles.get(p.id).setColor(Color.LIGHTSKYBLUE);

        }
        waitBetweenEmits++;
        if (particleCoordinates.size() < NUM_PARTICLES-1  && waitBetweenEmits > 100) {
            addParticles();
            waitBetweenEmits = 0;
            tree = new KdTree(2,particleCoordinates);
            tree.setRadius(smoothingLength);
            if (particleCoordinates.size() == NUM_PARTICLES){
                emitter = false;
            }
        }
    }
    private void addParticles(){
        Random rand = new Random();
        for (int i = 0; i < emmitParticlesNum; i++) {
            var p = new ParticleDrawn();
            particles.add(p);
            root.getChildren().add(p);
        }

        int x = emitterPosition[0];
        int y = emitterPosition[1];
        for (int i = particles.size()-emmitParticlesNum; i < particles.size(); i++) {
            ParticleDrawn particle = particles.get(i);
            double num = smoothingLength;

            particle.setTranslateX(x);
            particle.setTranslateY(y);
            positions[i][0] = x; //+Math.random()*1;
            positions[i][1] = y;



            x+= num/2;
        }

        for (int i = particles.size()-emmitParticlesNum; i < particles.size(); i++) {
            KdTree.Node newNode = new KdTree.Node(new double[] {positions[i][0],positions[i][1]});
            double [] velocity = {rand.nextInt(200)+1,rand.nextInt(200)+1};

            newNode.setVelocity(velocity);
            newNode.setDensity(BASE_DENSITY);
            newNode.id = i;
            particleCoordinates.add(newNode);

        }

    }

    private void onUpdateParalel(){



        List<KdTree.Node>[] points;
        points = new List[num_threads];
        for (int i = 0; i < points.length; i++) {
            points[i] = new ArrayList<>();
        }
            for (int i = 0;i < num_threads;i++) {
                for (int j = 0; j < particleCoordinates.size(); j++) {
                    if ((particleCoordinates.get(j).getCoords_()[0] > (MAX_X*i/num_threads) || (i==0 && particleCoordinates.get(j).getCoords_()[0] >= ((MAX_X*i/num_threads))) ) && particleCoordinates.get(j).getCoords_()[0] <= (MAX_X*(i+1)/num_threads)){
                        points[i].add(particleCoordinates.get(j));
                    }
                }
            }
            for (int i = 0; i < points.length; i++) {
                NeighbourThread nbt = new NeighbourThread(points[i]);
                executorService.submit(nbt);
            }

            executorService.shutdown();
            try {
                executorService.awaitTermination(Long.MAX_VALUE, TimeUnit.NANOSECONDS);
            } catch (InterruptedException e) {
                e.printStackTrace();
            }
            executorService = Executors.newCachedThreadPool();
            for (int i = 0; i < points.length; i++) {
                Segment segment = new Segment(points[i],tree);
                executorService.submit(segment);
            }

        /*for (int i = 0; i < points.length; i++) {
            for (int j = 0; j < points[i].size(); j++) {
                System.out.print(points[i].get(j).id+ ", ");
            }
            System.out.println();
        }

         */



        updateParticles();
    }

    private void onUpdate(){

        for (KdTree.Node particleCoordinate : particleCoordinates) {
            neighbors[particleCoordinate.id] = (tree.rangeSearch(particleCoordinate,smoothingLength));
            //neighbors[particleCoordinate.id] = findCollisionsINFLUID(particleCoordinate);
        }



        for (int i = 0; i < particles.size(); i++) {
            if (simulateDamBreak){
                if (dam < 1000){

                    damBreak(particleCoordinates.get(i));
                }
            }


            checkIfBounced(particleCoordinates.get(i));
        }
        if (simulateDamBreak){
            dam+=1;
        }

        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            calculateDensity(currParticle);
            if (currParticle.getDensity() == 0){
                currParticle.setDensity(0.0001);
            }
        }

        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            calculatePressure(currParticle);
        }

        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            calculateForces(currParticle);
        }

        // THIS IS THE SURFACE TENSION!

        List<KdTree.Node> surfaceParticles = getSurfaceParticles(particleCoordinates);

        List<double[]> normalVectors = getNormalVectors(surfaceParticles);

        List<Double> curvatureValues = findCurves(surfaceParticles, normalVectors);
        if (true){
            for (int i = 0;i < surfaceParticles.size();i++) {
                double[] force = calculateSurfaceTensionForce( curvatureValues.get(i),normalVectors.get(i),0.02, surfaceParticles.get(i).getPressure());
                if (add(surfaceParticles.get(i).getVelocity(),force)[1]+"" != "NaN") {
                    if (turnOnSurfaceTension){
                        newForces[i] = (add(surfaceParticles.get(i).getVelocity(),force));
                    }
                }
            }
        }
        for (int i = 0; i < particleCoordinates.size(); i++) {

            KdTree.Node currParticle = particleCoordinates.get(i);
            currParticle.setVelocity(newForces[i]);
            currParticle.setCoords_(add(currParticle.getCoords_(), mul(currParticle.getVelocity(),timeStep)));
        }

        // THIS IS THE SURFACE LINE CODE
        /*
        List<double[]> points = new ArrayList<>();
        for (int i = 0;i < surfaceParticles.size();i++) {
            points.add(surfaceParticles.get(i).getCoords_());
        }
        for (int i = 0; i < points.size(); i++) {
            for (int j = 0; j < points.size()-1; j++) {
                if (points.get(j)[0] > points.get(j+1)[0]){
                    Collections.swap(points,j,j+1);

                }
            }
        }
        List<double[]> surface_points = new ArrayList<>();
        for (int i = 0; i <= surface_detail; i++) {
            double[] largest = {MAX_X*i/surface_detail,0};
            for (int j = 0; j < points.size()-1; j++) {
                double[] curr = points.get(j);
                if (curr[0] > MAX_X*i/surface_detail && curr[0] > MAX_X*(i+1)/surface_detail){
                    if (largest[1] < curr[1]){
                        largest = curr;
                    }
                }

            }
            //System.out.println(largest[1] +" ---> " + MAX_X*(i+1)/surface_detail);
            surface_points.add(largest);
        }

         */
        //drawCurves(surface_points);

        updateParticles();

    }

    public void drawCurves(List<double[]> points){
        for (int i = 0; i < num_threads; i++) {
            Line curve = curves.get(i);
            //System.out.println(points.get(i)[0] + ", " +points.get(i)[1]);
            double[] Start = points.get(i);
            double[] End = points.get(i+1);


            curve.setStartX(Start[0]);
            curve.setStartY(MAX_Y-Start[1]);
            curve.setEndX(End[0]);
            curve.setEndY(MAX_Y-End[1]);



            curve.setStroke(Color.BROWN);
            curve.setStrokeWidth(3);

        }
        try {
            Thread.sleep(10);
        } catch (InterruptedException e) {
            e.printStackTrace();
        }
        System.out.println("_--------------------------------------");
    }


    public static double[] calculateSurfaceTensionForce(double curve,double[] normal, double surfaceTensionCoefficient, double pressure) {
        // Calculate the surface tension force using the Young-Laplace equation
        double[] force = new double[2];

        double normalLength = Math.sqrt(normal[0] * normal[0] + normal[1] * normal[1]);
        normal[0] = normal[0] / normalLength;
        normal[1] = normal[1] / normalLength;

        force[0] = surfaceTensionCoefficient * curve * normal[1];
        force[1] = -surfaceTensionCoefficient * curve * normal[0];

        return force;
    }


    public static List<Double> findCurves(List<KdTree.Node> surfaceParticles, List<double[]> NormalVectors){
        List<Double> curves = new ArrayList<>();
        for (int i = 0;i < surfaceParticles.size();i++) {
            double curvature = calculateCurvature(surfaceParticles.get(i), surfaceParticles, NormalVectors.get(i));

            // Add the curvature value to the list of curvature values
            curves.add(curvature);
        }
        return curves;
    }
    public static double calculateCurvature(KdTree.Node particle, List<KdTree.Node> surfaceParticles, double[] normal) {
        // Define a coordinate system centered at the particle
        double[] xAxis = { 1, 0 };
        double[] yAxis = { 0, 1 };
        double[] zAxis = normal;

        // Project the positions of the neighboring particles onto the xy-plane of the coordinate system
        List<double[]> projectedPositions = new ArrayList<>();
        for (KdTree.Node neighbor : surfaceParticles) {
            if (neighbor == particle) continue; // skip the current particle
            double[] projection = projectPosition(neighbor.getCoords_()[0], neighbor.getCoords_()[1], zAxis, xAxis, yAxis);
            projectedPositions.add(projection);
        }

        // Fit a curve to the projected positions using linear regression
        double[] coefficients = linearRegression(projectedPositions);

        // Calculate the curvature using the second derivative of the curve
        double curvature = 2 * coefficients[1];

        return curvature;
    }

    public static double[] projectPosition(double x, double y, double[] normal, double[] xAxis, double[] yAxis) {
        // Calculate the projection of the position onto the plane
        double[] projection = new double[2];
        projection[0] = x * xAxis[0] + y * yAxis[0];
        projection[1] = x * xAxis[1] + y * yAxis[1];

        return projection;
    }

    public static double[] linearRegression(List<double[]> points) {
        int n = points.size();

        // Calculate the sums
        double sx = 0, sy = 0, sxx = 0, syy = 0, sxy = 0;
        for (double[] point : points) {
            double x = point[0];
            double y = point[1];
            sx += x;
            sy += y;
            sxx += x * x;
            syy += y * y;
            sxy += x * y;
        }

        // Calculate the coefficients
        double a = (n * sxy - sx * sy) / (n * sxx - sx * sx);
        double b = (sy - a * sx) / n;

        return new double[] { a, b };
    }

    public static List<double[]> getNormalVectors(List<KdTree.Node> surfaceParticles){
        List<double[]> normalVectors = new ArrayList<>();
        for ( KdTree.Node particle : surfaceParticles){
            double[] avgNormal = new double[2];
            for (KdTree.Node neighbour : neighbors[particle.id]){
                double[] faceNormal = calculateFaceNormal(particle.getCoords_()[0],particle.getCoords_()[1],neighbour.getCoords_()[0],neighbour.getCoords_()[1]);
                avgNormal[0] += faceNormal[0];
                avgNormal[1] += faceNormal[1];
            }
            double length = Math.sqrt(avgNormal[0] * avgNormal[0] + avgNormal[1] * avgNormal[1]);
            avgNormal[0] /= length;
            avgNormal[1] /= length;

            normalVectors.add(avgNormal);
        }
        return normalVectors;
    }
    public static double[] calculateFaceNormal(double x1, double y1, double x2, double y2) {
        // Calculate the vector pointing from (x1, y1) to (x2, y2)
        double[] v = { x2 - x1, y2 - y1 };

        // Rotate the vector 90 degrees clockwise to obtain the normal vector
        double[] n = { -v[1], v[0] };

        return n;
    }
    public static List<KdTree.Node> getSurfaceParticles(List<KdTree.Node> findArray){
        List<KdTree.Node> surfaceParticles = new ArrayList<>();
        for (int i = 0; i < findArray.size(); i++) {
            KdTree.Node currParticle = findArray.get(i);
            if (find_surface(currParticle)){
                surfaceParticles.add(currParticle);
            }
        }
        return surfaceParticles;
    }

    public static boolean find_surface(KdTree.Node currParticle){
        double densityThreshold = -19.94;
        if (PARTICLE_RADIUS < 6){
            densityThreshold = -19.5;
        }
        //System.out.println(currParticle.getPressure());

        if (currParticle.getPressure() < densityThreshold&& neighbors[currParticle.id].size() < 3) {
            if (colorfull){

            }
            currParticle.setColor(Color.WHITE);
            if (decrease_foam_volume){
                currParticle.setRadius(PARTICLE_RADIUS*2-2);
            }
            if(neighbors[currParticle.id].size() > 1){
                return true;
            }
        }else {
            if (neighbors[currParticle.id].size() < 5){
                if (colorfull){
                    currParticle.setColor(Color.LIGHTSKYBLUE);
                }
                //
                if (decrease_foam_volume){
                    currParticle.setRadius(PARTICLE_RADIUS*2);
                }
            }
            else {
                //
                if (colorfull){
                    currParticle.setColor(Color.DEEPSKYBLUE);
                }
                if (decrease_foam_volume){
                    currParticle.setRadius(PARTICLE_RADIUS*2+2);
                }
            }
        }
        return false;
    }

    public static void main(String[] args) {

        executorService = Executors.newWorkStealingPool();

        launch();



    }


    public static double getDistance(double particleX, double particleY, double neighbourX, double neighbourY){
        double dx = particleX - neighbourX;
        double dy = particleY - neighbourY;
        double distance = Math.sqrt(dx * dx + dy * dy);

        if (distance <= 0){
            distance = 0.001;
        }
        return distance;
    }


}
