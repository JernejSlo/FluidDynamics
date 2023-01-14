package com.example.fluiddynamicsjavafx;

import javafx.animation.AnimationTimer;
import javafx.application.Application;
import javafx.scene.Parent;
import javafx.scene.Scene;
import javafx.scene.layout.Pane;
import javafx.scene.paint.Color;
import javafx.scene.shape.Circle;
import javafx.stage.Stage;

import java.io.IOException;
import java.util.*;


public class Fluid2D extends Application {

    int SCREEN_WIDTH = 1280;
    int SCREEN_HEIGHT = 720;
    int dam = 0;
    boolean simulateDamBreak = false;


    public static KdTree tree;
    public static int NUM_PARTICLES = 2000;
    ArrayList<KdTree.Node> [] neighbors = new ArrayList[NUM_PARTICLES];
    public static int PARTICLE_RADIUS = 8;
    private List<ParticleDrawn> particles = new ArrayList<>();

    public static double [][] positions = new double[NUM_PARTICLES][2];
    public static double [] densities = new double[NUM_PARTICLES];

    public static List<KdTree.Node> particleCoordinates = new ArrayList<>();

    double smoothingLength = PARTICLE_RADIUS*2;



    double mass = 1;
    double [] zeros;
    double timeStep = 0.01;

    public static double ISOTROPIC_EXPONENT = 20;
    public static double BASE_DENSITY = 1;
    public static double [] CONSTANT_FORCE = {0,-0.1};
    public static double DYNAMIC_VISCOSITY = 0.5;
    public static double DAMPING_COEFFICIENT = - 0.9;

    public double NORMALIZATION_DENSITY = (315 * mass) / (64 * Math.PI * Math.pow((smoothingLength),9));
    public double NORMALIZATION_VISCOUS_FORCE = (45 * DYNAMIC_VISCOSITY * mass) / (Math.PI * (Math.pow((smoothingLength),6)));
    public double NORMALIZATION_PRESSURE_FORCE = -((45 * mass) / (Math.PI * Math.pow((smoothingLength),6)));

    // Define the time step for the simulation
    double dt = 0.001;

    private static final int MIN_X = 0;
    private static final int MIN_Y = 0;
    private double MAX_X = SCREEN_WIDTH;
    private double MAX_Y = SCREEN_HEIGHT;

    // The gravitational force acting on the vectors
    private static final double GRAVITY = 9.8;

    // The coefficient of restitution for the vectors
    private static final double COR = 0.9;

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
            double [] velocity = {400,-200};

            newNode.setVelocity(velocity);
            newNode.setDensity(BASE_DENSITY);
            newNode.id = i;
            particleCoordinates.add(newNode);

        }

        tree = new KdTree(2,particleCoordinates);
        tree.setRadius(smoothingLength);



        AnimationTimer timer = new AnimationTimer() {
            @Override
            public void handle(long now) {
                onUpdate();

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



    public void calculateDensity(KdTree.Node currParticle){
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

    public void calculatePressure(KdTree.Node currParticle){
        double pressure = ISOTROPIC_EXPONENT * (currParticle.getDensity() -BASE_DENSITY);
        currParticle.setPressure(pressure);
    }

    public double[] normalize(double[] vector){
        double magnitude = Math.sqrt(vector[0]*vector[0] + vector[1]*vector[1]);
        double[] normalizedVector = {vector[0]/magnitude, vector[1]/magnitude};



        return normalizedVector;
    }
    public double [] sub(double [] vector,double[] subtract){
        double [] subbedVector = new double[2];

        subbedVector[0] = vector[0]-subtract[0];
        subbedVector[1] = vector[1]-subtract[1];

        return subbedVector;
    }
    public double [] add(double [] vector,double[] add){
        double [] addedVector = new double[2];

        addedVector[0] = vector[0]+add[0];
        addedVector[1] = vector[1]+add[1];

        return addedVector;
    }
    public double [] mul(double [] vector1,double scalar){
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
    public double [] div(double [] vector1,double scalar){
        double [] divVector = new double[2];

        divVector[0] = vector1[0]/scalar;
        divVector[1] = vector1[1]/scalar;

        return divVector;
    }
    public double dot(double [] vector1,double [] vector2){

        double dotProduct = (vector1[0]* vector2[0]) + (vector1[1]* vector2[1]);

        return dotProduct;
    }

    public void calculateForces(KdTree.Node currParticle) {
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
        currParticle.setVelocity(velocityCalculated);
        currParticle.setCoords_(add(pl,mul(currParticle.getVelocity(),timeStep)));

    }

    public void damBreak(KdTree.Node currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] >= MAX_X/2+smoothingLength){

            double [] newVelocity = {-velocity[0]*COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {MAX_X/2+smoothingLength, currparticle.getCoords_()[1]});
        }
    }

    // dodaj da se odbije od zida
    public void checkIfBounced(KdTree.Node currparticle){
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

        }
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

        List<KdTree.Node> surfaceParticles = getSurfaceParticles();
        /*
        List<double[]> normalVectors = getNormalVectors(surfaceParticles);

        List<Double> curvatureValues = findCurves(surfaceParticles, normalVectors);

        for (int i = 0;i < surfaceParticles.size();i++) {
            double[] force = calculateSurfaceTensionForce( curvatureValues.get(i),normalVectors.get(i),0.01, surfaceParticles.get(i).getPressure());
            if (add(surfaceParticles.get(i).getVelocity(),force)[1]+"" != "NaN") {
                surfaceParticles.get(i).setVelocity(add(surfaceParticles.get(i).getVelocity(),force));
            }
        }

         */





        updateParticles();

    }

    public static double[] calculateSurfaceTensionForce(double curve,double[] normal, double surfaceTensionCoefficient, double pressure) {
        // Calculate the surface tension force using the Young-Laplace equation
        double[] force = new double[2];
        force[0] = surfaceTensionCoefficient * curve * normal[1];
        force[1] = -surfaceTensionCoefficient * curve * normal[0];

        return force;
    }


    public List<Double> findCurves(List<KdTree.Node> surfaceParticles,List<double[]> NormalVectors){
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

    public List<double[]> getNormalVectors(List<KdTree.Node> surfaceParticles){
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
    public List<KdTree.Node> getSurfaceParticles(){
        List<KdTree.Node> surfaceParticles = new ArrayList<>();
        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            if (find_surface(currParticle)){
                surfaceParticles.add(currParticle);
            }
        }
        return surfaceParticles;
    }

    public boolean find_surface(KdTree.Node currParticle){
        double densityThreshold = -19.94;
        System.out.println(currParticle.getPressure());

        if (currParticle.getPressure() < densityThreshold&& neighbors[currParticle.id].size() < 3) {
            currParticle.setColor(Color.WHITE);
            currParticle.setRadius(PARTICLE_RADIUS*2-1);
            return true;
        }else {
            if (neighbors[currParticle.id].size() < 5){
                currParticle.setColor(Color.LIGHTSKYBLUE);
                currParticle.setRadius(PARTICLE_RADIUS*2);
            }
            else {
                currParticle.setColor(Color.DEEPSKYBLUE);
                currParticle.setRadius(PARTICLE_RADIUS*2+2);
            }
        }
        return false;
    }

    public static void main(String[] args) {
        launch();
    }


    public double getDistance(double particleX,double particleY,double neighbourX,double neighbourY){
        double dx = particleX - neighbourX;
        double dy = particleY - neighbourY;
        double distance = Math.sqrt(dx * dx + dy * dy);

        if (distance <= 0){
            distance = 0.001;
        }
        return distance;
    }


}