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


public class HelloApplication extends Application {

    int SCREEN_WIDTH = 720;
    int SCREEN_HEIGHT = 720;
    public static int NUM_PARTICLES = 100;
    public static int PARTICLE_RADIUS = 1;
    private List<ParticleDrawn> particles = new ArrayList<>();

    public static double [][] positions = new double[NUM_PARTICLES][2];
    public static double [] densities = new double[NUM_PARTICLES];

    public static List<KdTree.Node> particleCoordinates = new ArrayList<>();
    public static KdTree tree;




    // ------------------------------------------------- //

    double smoothingLength = 2*PARTICLE_RADIUS;
    double mass = 0.15;
    double stiffness = 100;
    double viscosity = 0.0001;
    double restDensity = 0.001; // rest density of the fluid
    double gasConstant = 0.0000001; // gas constant of the fluid
    double [] zeros;
    double timeStep = 8;

    // Define the time step for the simulation
    double dt = 0.001;

    private static final int MIN_X = 0;
    private static final int MIN_Y = 0;
    private double MAX_X = SCREEN_WIDTH;
    private double MAX_Y = SCREEN_HEIGHT;

    // The gravitational force acting on the vectors
    private static final double GRAVITY = 9.8;

    // The coefficient of restitution for the vectors
    private static final double COR = 0.8;

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


        int x = 10;
        int y = 10;
        for (int i = 0; i < NUM_PARTICLES; i++) {
            ParticleDrawn particle = particles.get(i);
            int num = 4;

            particle.setTranslateX(x);
            particle.setTranslateY(y);
            positions[i][0] = x;
            positions[i][1] = y;

            if (x >= MAX_X-10){
                x = 0;
                y += num;
            }
            x+=num;

                /*
                particle.setTranslateX(10);
                particle.setTranslateY(10);

                positions[i*50+j][0] = 10;
                positions[i*50+j][1] = 10;

                 */
        }

        for (int i = 0; i < positions.length; i++) {
            KdTree.Node newNode = new KdTree.Node(new double[] {positions[i][0],positions[i][1]});
            double [] velocity = {1,-1};

            newNode.setVelocity(velocity);
            newNode.setDensity(0.1);
            newNode.id = i;
            particleCoordinates.add(newNode);
            newNode.show();

        }

        tree = new KdTree(2,particleCoordinates);


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
        ParticleDrawn(){
            var shape = new Circle(PARTICLE_RADIUS*2);
            shape.setFill(Color.LIGHTSKYBLUE);

            getChildren().addAll(shape);
        }
    }


    double calculateKernel(double r,double smoothingLength) {
        if (r <= 0) return 0;
        double sigma = smoothingLength; // kernel radius
        double q = r / sigma;
        if (q <= 1) {
            return 15.0 / (7.0 * Math.PI * sigma * sigma) * (1.0 - 1.5 * q * q + 0.75 * q * q * q);
        } else if (q <= 2) {
            return 15.0 / (7.0 * Math.PI * sigma * sigma) * (0.25 * (2.0 - q) * (2.0 - q) * (2.0 - q));
        } else {
            return 0.1;
        }
    }

    public void calculateDensity(KdTree.Node currParticle){
        ArrayList collisions = tree.findCollisions(currParticle, smoothingLength);
        double [] pl = currParticle.getCoords_();
        double density = 0;
        //System.out.println("before density");
        //currParticle.show();
        for (int i = 0; i < collisions.size(); i++) {
            KdTree.Node neighbour = (KdTree.Node) collisions.get(i);
            double [] c = neighbour.getCoords_();
            double distance = getDistance(pl[0],pl[1],c[0],c[1]);
            double kernelValue = calculateKernel(distance,smoothingLength);
            //System.out.println("kernel " + kernelValue);
            //System.out.println("particle density " + currParticle.getDensity());
            //System.out.println("smoothingLen " + smoothingLength);
            density+=mass*kernelValue;
            //System.out.println("step " + i +" -> density: "+density);

        }
        currParticle.setDensity(density);
        //System.out.println("final density" + density);
        //System.out.println("after density");
        //currParticle.show();



    }

    public void calculatePressure(KdTree.Node currParticle){
        double pressure = gasConstant * (currParticle.getDensity() - restDensity);
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

    public double dot(double [] vector1,double [] vector2){

        double dotProduct = (vector1[0]* vector2[0]) + (vector1[1]* vector2[1]);

        return dotProduct;
    }

    public void calculateForces(KdTree.Node currParticle) {
        double[] velocity = currParticle.getVelocity();
        double[] pressureVelocity = {0,0};
        double[] viscousVelocity = {0,0};
        /*
        System.out.println("before forces");
        currParticle.show();
        System.out.println("");

         */
        ArrayList collisions = tree.findCollisions(currParticle, smoothingLength);
        double [] pl = currParticle.getCoords_();
        for (int i = 0; i < collisions.size(); i++) {
            KdTree.Node neighbour = (KdTree.Node) collisions.get(i);
            double [] c = neighbour.getCoords_();
            double distance = getDistance(pl[0],pl[1],c[0],c[1]);
            if (distance == 0){
                distance = 0.01;
            }
            double kernelValue = calculateKernel(distance,smoothingLength);
            double P = currParticle.getPressure() + neighbour.getPressure();
            double mu = 2*viscosity;

            double [] rhat = normalize(sub(currParticle.getCoords_(),neighbour.getCoords_()));

            pressureVelocity = add(pressureVelocity,mul(rhat,(P*kernelValue)));

            double[] v = sub(currParticle.getVelocity(), neighbour.getVelocity());
            viscousVelocity = add(viscousVelocity,mul(rhat,(-mu*(dot(rhat,v)*kernelValue))));

            /*
            System.out.println("v: " + v[0]+", "+v[1]);
            System.out.println("rhat: " + rhat[0]+", "+rhat[1]);
            System.out.println("mu: " + mu);
            System.out.println("P: " + P);
            System.out.println("kernelValue: " + kernelValue);
            System.out.println("pressureVelocity: " + pressureVelocity[0]+", "+pressureVelocity[1]);
            System.out.println("addition to pressureVelocity: " + mul(rhat,(P*kernelValue))[0]+", "+mul(rhat,(P*kernelValue))[1]);
            System.out.println("-----------------------------------------------------");

             */






        }
        //System.out.println("velocity: "+ velocity[0] + " " + velocity[1]);
        double gforce = 0.0000098;
        double [] force = sub(add(pressureVelocity,viscousVelocity), new double[]{0, gforce*timeStep});

        System.out.println("Force: "+force[0] + " " + force[1]);

        currParticle.setVelocity(add(currParticle.getVelocity(),mul(force,(timeStep/mass))));
        currParticle.setCoords_(add(currParticle.getCoords_(), mul(currParticle.getVelocity(), timeStep)));

        /*
        System.out.println("after forces");
        currParticle.show();
        System.out.println();

         */
    }

    // dodaj da se odbije od zida
    public void checkIfBounced(KdTree.Node currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] <= MIN_X){

            System.out.println("bounced off of left wall");
            double [] newVelocity = {-velocity[0]*COR,velocity[1]*COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {0, currparticle.getCoords_()[1]});
        }
        if (coords[0] >= MAX_X){
            System.out.println("bounced off of right wall");

            double [] newVelocity = {-velocity[0]*COR,velocity[1]*COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {MAX_X, currparticle.getCoords_()[1]});
        }
        if (coords[1] <= 0){
            System.out.println("bounced off of floor");
            //double [] newVelocity = {velocity[0],-velocity[1]};
            double [] newVelocity = {velocity[0]*COR,-velocity[1]*COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],0});
        }
        if (coords[1] >= MAX_Y){

            System.out.println("bounced off of celling");
            double [] newVelocity = {velocity[0]*COR,-velocity[1]*COR};
            //double [] newVelocity = {velocity[0],-velocity[1]};
            currparticle.setVelocity(newVelocity);

            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],MAX_Y});
        }
    }

    public void updateParticles() {

        for (KdTree.Node p : particleCoordinates) {


            particles.get(p.id).setTranslateX(p.getCoords_()[0]);
            particles.get(p.id).setTranslateY(Math.abs(p.getCoords_()[1]-MAX_Y));
            p.show();

            /*
            particles.get(p.id).setTranslateX(p.getCoords_()[0]);
            particles.get(p.id).setTranslateY(p.getCoords_()[1]);

             */

        }
    }

    private void onUpdate(){

        KdTree.Node nearest;

        for (int i = 0; i < particles.size(); i++) {
            checkIfBounced(particleCoordinates.get(i));
        }

        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            calculateDensity(currParticle);
            if (currParticle.getDensity() == 0){
                currParticle.setDensity(0.0001);
            }
            //System.out.println("in density calc "+currParticle.getVelocity()[0]+" "+currParticle.getVelocity()[1]);
        }

        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            calculatePressure(currParticle);
        }

        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            calculateForces(currParticle);
            //System.out.println("in forces calc "+currParticle.getVelocity()[0]+" "+currParticle.getVelocity()[1]);
        }


        updateParticles();


        tree = new KdTree(2,particleCoordinates);

    }

    public static void main(String[] args) {
        launch();
    }


    public double getDistance(double particleX,double particleY,double neighbourX,double neighbourY){
        double distance = 0;
        double dx = Math.abs(neighbourX - particleX);
        double dy = Math.abs(neighbourY - particleY);

        double min = Math.min(dx, dy);
        double max = Math.max(dx, dy);

        double diagonalSteps = min;
        double straightSteps = max - min;

        distance = Math.sqrt(2) * diagonalSteps + straightSteps;


        return distance;
    }


}
