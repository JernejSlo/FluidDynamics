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
import java.util.concurrent.CyclicBarrier;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;


public class Fluid2D extends Application {

    // Project by : Jernej Koprivnikar
    //
    // Credit to the creator of this video - https://www.youtube.com/watch?v=-0m05gzk8nk&ab_channel=MachineLearning%26Simulation
    // some equations used in my code are from this video (density, pressure and forces. The original code is in python).
    // The creator of the video used the calculations from this research paper to create his code:
    // https://matthias-research.github.io/pages/publications/sca03.pdf


    // CHANGEABLE PARAMETERS

    public static int SCREEN_WIDTH = 1200;
    public static int SCREEN_HEIGHT = 800;

    // enable distributed mode
    // to simulate start the server with Fluid2D and a client with test and client
    public boolean distributed = false;
    public boolean host = true;
    public static String hostAddress = "127.0.0.1";
    public static int hostPort = 9999;

    public static double prev = System.nanoTime();
    public static double[] all = {0,0};

    // run it in parallel mode
    public boolean parallel = true;

    // turns emitter mode on/off
    public boolean emitter = false;
    // location of the emitter, don't place it outside the window!
    public int[] emitterPosition = {SCREEN_WIDTH/2,SCREEN_HEIGHT/2};
    // time between emits
    public int waitBetweenEmits = 10;

    // when timeStep is too large it moves the particles too fast and creates a buggy look
    // it performs the best with slower time steps      !This is true for all modes, since the calculations are more precise!
    public static boolean turnOnSurfaceTension = false;

    // adds color to the particles
    public static boolean colorful = true;

    // running with this turned off is somewhat ugly
    // this increases the size of all particles by a varying amount based on the pressure
    public static boolean decrease_foam_volume = true;

    // increasing this means shortening the time that the particles are in the "dam"
    static int dam = 0;
    // this starts the simulation by simulating a dam break scenario
    static boolean simulateDamBreak = false;

    // time step can be changed, with lower values (0.005>) it works better, but it is very slow
    // also SLOWER timeStep show the surface tension way better, you can clearly see how particles join together to form a surface
    // values that are greater than 0.01, can sometimes mess up the calculations but most of the time under 0.02 is fine
    static double timeStep = 0.015;


    public static int PARTICLE_RADIUS = 8;
    public static int NUM_PARTICLES = 2000;
    protected static ArrayList<Particle> [] neighbors = new ArrayList[NUM_PARTICLES];
    double[] startVector = {500,-200};


    // SOME CAN BE CHANGED BUT REALLY SHOULDN'T BE
    public static ExecutorService executorService;

    public int emmitParticlesNum = 10;
    public static Grid grid;
    public static int size = 12;
    public static int num_threads = 48;

    private static List<ParticleDrawn> particles = new ArrayList<>();
    private List<Line> curves = new ArrayList<>();
    public static double[][] newForces = new double[NUM_PARTICLES][2];

    public static double [][] positions = new double[NUM_PARTICLES][2];
    public static List<Particle> particleCoordinates = new ArrayList<>();

    static double smoothingLength = PARTICLE_RADIUS*2;
    Pane root;

    public int emmitCounter = 10;
    static double mass = 1;
    double [] zeros;

    public static double ISOTROPIC_EXPONENT = 20;
    public static double BASE_DENSITY = 1;
    public static double [] CONSTANT_FORCE = {0.05,-0.1};
    public static double DYNAMIC_VISCOSITY = 0.5;

    public static double NORMALIZATION_DENSITY = (315 * mass) / (64 * Math.PI * Math.pow((smoothingLength),9));
    public static double NORMALIZATION_VISCOUS_FORCE = (45) / (Math.PI * (Math.pow((smoothingLength),6)));
    public static double NORMALIZATION_PRESSURE_FORCE = (-(45) / (Math.PI * Math.pow((smoothingLength),6)));

    public static final int MIN_X = 0;
    private static final int MIN_Y = 0;
    public static double MAX_X = SCREEN_WIDTH;
    protected static double MAX_Y = SCREEN_HEIGHT;

    // The coefficient of restitution for the vectors
    public static final double COR = 0.9;
    private boolean started = false;

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
                    positions[i][0] = x;
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
                    positions[i][0] = x;
                    positions[i][1] = y;

                    if (x >= (3*MAX_X)/4){
                        x = MAX_X/4;
                        y += num;
                    }

                    x+= num/2;

                }
            }
            for (int i = 0; i < positions.length; i++) {
                Particle newNode = new Particle(new double[] {positions[i][0],positions[i][1]},i);
                double [] velocity = startVector;

                newNode.setVelocity(velocity);
                newNode.setDensity(BASE_DENSITY);
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
                positions[i][0] = x;
                positions[i][1] = y;



                x+= num/2;

            }

            for (int i = 0; i < emmitParticlesNum; i++) {
                Particle newNode = new Particle(new double[] {positions[i][0],positions[i][1]},i);
                double [] velocity = {rand.nextInt(200)+1,rand.nextInt(200)+1};

                newNode.setVelocity(velocity);
                newNode.setDensity(BASE_DENSITY);
                particleCoordinates.add(newNode);

            }
        }




        grid = new Grid(size,smoothingLength);
        for (Particle p : particleCoordinates){
            grid.Insert(p);
        }



        AnimationTimer timer = new AnimationTimer() {
            @Override
            public void handle(long now) {
                if (distributed){
                    onUpdateDistributed();
                }
                else {
                    if (parallel){
                        onUpdateParallel();
                    }
                    else {

                        onUpdate();
                    }
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

    public static void calculateDensity(Particle currParticle){
        ArrayList collisions = neighbors[currParticle.id];
        double [] pl = currParticle.getCoords_();
        double density = Math.pow((Math.pow(smoothingLength,2)-Math.pow(0,2)),3);
        for (int i = 0; i < collisions.size(); i++) {
            Particle neighbour = (Particle) collisions.get(i);
            double [] c = neighbour.getCoords_();
            double distance = getDistance(pl[0],pl[1],c[0],c[1]);
            density+=Math.pow((Math.pow(smoothingLength,2)-Math.pow(distance-PARTICLE_RADIUS,2)),3);

        }
        currParticle.setDensity(NORMALIZATION_DENSITY * density);
    }


    public static void calculatePressure(Particle currParticle){
        double pressure = ISOTROPIC_EXPONENT * (currParticle.getDensity() -BASE_DENSITY);
        currParticle.setPressure(pressure);
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
    public static double [] div(double[] vector1, double scalar){
        double [] divVector = new double[2];

        divVector[0] = vector1[0]/scalar;
        divVector[1] = vector1[1]/scalar;

        return divVector;
    }

    public static void calculateForces(Particle currParticle) {
        double[] velocity = {0,0};
        double[] pressureForce = {0,0};
        double[] viscousForce = {0,0};
        double[] currentVelocity = currParticle.getVelocity();
        ArrayList collisions = neighbors[currParticle.id];
        double [] pl = currParticle.getCoords_();
        for (int i = 0; i < collisions.size(); i++) {
            Particle neighbour = (Particle) collisions.get(i);
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

    }

    public static void damBreak(Particle currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] >= MAX_X/2+smoothingLength){

            double [] newVelocity = {-velocity[0]*COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {MAX_X/2+smoothingLength, currparticle.getCoords_()[1]});
        }
    }

    // check if the particle has bounced from the walls and adjust vector speed
    public static void checkIfBounced(Particle currparticle){
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
            double [] newVelocity = {velocity[0],-velocity[1]*COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],0});
        }
        if (coords[1] >= MAX_Y){

            double [] newVelocity = {velocity[0],-velocity[1]*COR};
            currparticle.setVelocity(newVelocity);

            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],MAX_Y});
        }
    }

    public static void updateParticlesDistributed() {

        double start1 = System.nanoTime();
        double last = (start1-prev)/1000000.0;
        System.out.println("Time since last update " + last);
        prev = start1;
        all[0]+= last;
        all[1]++;
        double calc = (all[0]/all[1]);
        System.out.println("Avg " + calc);

        for (Particle p : particleCoordinates) {
            particles.get(p.id).setTranslateX(p.coords_[0]);
            particles.get(p.id).setTranslateY(Math.abs(p.getCoords_()[1]-MAX_Y));
            particles.get(p.id).setColor(p.getColor());
            particles.get(p.id).setRadius(p.getRadius());
        }

        grid = new Grid(size,smoothingLength);
        for (Particle p : particleCoordinates){
            grid.Insert(p);
        }

    }
    // updates the positions of the particles on screen and takes care of adding new particles for emmiter
    public void updateParticles() {

        double start1 = System.nanoTime();
        double last = (start1-prev)/1000000.0;
        System.out.println("Time since last update " + last);
        prev = start1;
        all[0]+= last;
        all[1]++;
        double calc = (all[0]/all[1]);
        System.out.println("Avg " + calc);


        for (Particle p : particleCoordinates) {



            particles.get(p.id).setTranslateX(p.coords_[0]);
            particles.get(p.id).setTranslateY(Math.abs(p.getCoords_()[1]-MAX_Y));
            particles.get(p.id).setColor(p.getColor());
            particles.get(p.id).setRadius(p.getRadius());


        }

        emmitCounter++;
        if (particleCoordinates.size() < NUM_PARTICLES-emmitParticlesNum  && emmitCounter > waitBetweenEmits) {
            addParticles();
            emmitCounter = 0;
            grid = new Grid(size,smoothingLength);
            for (Particle p : particleCoordinates){
                grid.Insert(p);
            }
            if (particleCoordinates.size() == NUM_PARTICLES){
                emitter = false;
            }
        }
        else {
            grid = new Grid(size,smoothingLength);
            for (Particle p : particleCoordinates){
                grid.Insert(p);
            }
        }
    }
    // adds the particle to the panel and other arrays
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
            positions[i][0] = x;
            positions[i][1] = y;



            x+= num/2;
        }

        for (int i = particles.size()-emmitParticlesNum; i < particles.size(); i++) {
            Particle newNode = new Particle(new double[] {positions[i][0],positions[i][1]},i);
            double [] velocity = {rand.nextInt(200)+1,rand.nextInt(200)+1};

            newNode.setVelocity(velocity);
            newNode.setDensity(BASE_DENSITY);
            particleCoordinates.add(newNode);

        }

    }

    // takes care of connecting, managing connections and everything to do with the distributed approach
    private void onUpdateDistributed(){

        if (!started){
            if (host){
                Thread serverThread = new Thread(() -> {
                    Server server = new Server(particleCoordinates);
                    server.run();
                });
                serverThread.start();
                started = true;
            }
            else {
                new Thread(new Client()).start();
            }
        }

    }

    // takes care of updating the particles for the parallel version
    private void onUpdateParallel(){

        List<Particle>[] points;
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

        CyclicBarrier barrier = new CyclicBarrier(num_threads);

        for (int i = 0; i < points.length; i++) {
            Segment segment = new Segment(points[i],grid,barrier);
            executorService.submit(segment);
        }

        updateParticles();
    }

    // main function that calls other functions for calculations
    private void onUpdate(){
        int size_ = particles.size();
        // actual function
        for (Particle particleCoordinate : particleCoordinates) {
            neighbors[particleCoordinate.id] = grid.GetNeighbors(particleCoordinate);

        }

        for (int i = 0; i < size_; i++) {
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

        for (int i = 0; i < size_; i++) {
            Particle currParticle = particleCoordinates.get(i);
            calculateDensity(currParticle);
            if (currParticle.getDensity() == 0){
                currParticle.setDensity(0.0001);
            }
        }

        for (int i = 0; i < size_; i++) {
            Particle currParticle = particleCoordinates.get(i);
            calculatePressure(currParticle);
        }

        for (int i = 0; i < size_; i++) {
            Particle currParticle = particleCoordinates.get(i);
            calculateForces(currParticle);
        }

        // THIS IS THE SURFACE TENSION!

        List<Particle> surfaceParticles = getSurfaceParticles(particleCoordinates);

        List<double[]> normalVectors = getNormalVectors(surfaceParticles);

        List<Double> curvatureValues = findCurves(surfaceParticles, normalVectors);
        if (true){
            for (int i = 0;i < surfaceParticles.size();i++) {
                double[] force = calculateSurfaceTensionForce( curvatureValues.get(i),normalVectors.get(i),0.0000000001, surfaceParticles.get(i).getPressure());
                if (add(surfaceParticles.get(i).getVelocity(),force)[1]+"" != "NaN") {
                    if (turnOnSurfaceTension){
                        newForces[i] = (add(newForces[surfaceParticles.get(i).id], force));
                    }
                }
            }
        }
        for (int i = 0; i < particleCoordinates.size(); i++) {

            Particle currParticle = particleCoordinates.get(i);
            currParticle.setVelocity(newForces[i]);
            currParticle.setCoords_(add(currParticle.getCoords_(), mul(currParticle.getVelocity(),timeStep)));
        }



        updateParticles();

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


    public static List<Double> findCurves(List<Particle> surfaceParticles, List<double[]> NormalVectors){
        List<Double> curves = new ArrayList<>();
        for (int i = 0;i < surfaceParticles.size();i++) {
            double curvature = calculateCurvature(surfaceParticles.get(i), surfaceParticles, NormalVectors.get(i));

            // Add the curvature value to the list of curvature values
            curves.add(curvature);
        }
        return curves;
    }
    public static double calculateCurvature(Particle particle, List<Particle> surfaceParticles, double[] normal) {
        // Define a coordinate system centered at the particle
        double[] xAxis = { 1, 0 };
        double[] yAxis = { 0, 1 };
        double[] zAxis = normal;

        // Project the positions of the neighboring particles onto the xy-plane of the coordinate system
        List<double[]> projectedPositions = new ArrayList<>();
        for (Particle neighbor : surfaceParticles) {
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

    public static List<double[]> getNormalVectors(List<Particle> surfaceParticles){
        List<double[]> normalVectors = new ArrayList<>();
        for ( Particle particle : surfaceParticles){
            double[] avgNormal = new double[2];
            for (Particle neighbour : neighbors[particle.id]){
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
    public static List<Particle> getSurfaceParticles(List<Particle> findArray){
        List<Particle> surfaceParticles = new ArrayList<>();
        for (int i = 0; i < findArray.size(); i++) {
            Particle currParticle = findArray.get(i);
            if (find_surface(currParticle)){
                surfaceParticles.add(currParticle);
            }
        }
        return surfaceParticles;
    }
    // finds all particles that are on the surface, judging by the number of its neighbours and its pressure. It also resizes and colors them.
    public static boolean find_surface(Particle currParticle){
        double densityThreshold = -19.94;
        if (PARTICLE_RADIUS < 6){
            densityThreshold = -19.5;
        }

        if (currParticle.getPressure() < densityThreshold&& neighbors[currParticle.id].size() < 3) {
            if (colorful){

            }
            currParticle.setColor(Color.WHITE);
            if (decrease_foam_volume){
                currParticle.setRadius(PARTICLE_RADIUS*2-2);
            }
            if(neighbors[currParticle.id].size() > 0){
                return true;
            }
        }else {
            if (neighbors[currParticle.id].size() < 5){
                if (colorful){
                    currParticle.setColor(Color.LIGHTSKYBLUE);
                }
                if (decrease_foam_volume){
                    currParticle.setRadius(PARTICLE_RADIUS*2);
                }
            }
            else {
                if (colorful){
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

    // distance function
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
