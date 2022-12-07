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

    int SCREEN_WIDTH = 400;
    int SCREEN_HEIGHT = 720;
    public static int NUM_PARTICLES = 20;
    public static int PARTICLE_RADIUS = 1;
    private List<Particle> particles = new ArrayList<>();
    public double [] vectorsX;
    public static double GRAVITATIONAL_FORCE = 0.02;
    public static double PARTICLE_WEIGHT = 0.000018;
    public static double MOMENTUM_LOSS =   0.000000005;
    public double [] vectorsY;

    //tutorial

    public static double MAX_PARTICLES = 125;
    public static double DOMAIN_WIDTH = 40;
    public static double DOMAIN_HEIGHT = 80;

    public static double PARTICLE_MASS = 1;
    public static double ISOTROPIC_EXPONENT = 20;
    public static double BASE_DENSITY = 1;
    public static double SMOOTHING_LENGTH = 5;
    public static double DYNAMIC_VISCOSITY = 0.5;
    public static double DAMPING_COEFFICIENT = - 0.9;
    public static double NORMALIZATION_DENSITY = (315 * PARTICLE_MASS) / (Math.pow((64 * Math.PI * SMOOTHING_LENGTH),9));
    public static double NORMALIZATION_PRESSURE_FORCE = (45 * DYNAMIC_VISCOSITY * PARTICLE_MASS) / (Math.pow((Math.PI * SMOOTHING_LENGTH),6));
    public static double NORMALIZATION_VISCOUS_FORCE = (45 * PARTICLE_MASS) / (Math.pow((Math.PI * SMOOTHING_LENGTH),6));

    public static double [] CONSTANT_FORCE = new double[2];

    public static double TIME_STEP_LENGTH = 0.01;
    public static double N_TIME_STEPS = 2500;
    public static double ADD_PARTICLES_EVERY = 0;

    public static double [] FIGURE_SIZE = new double[2];
    public static double PLOT_EVERY = 6;
    public static double SCATTER_DOT_SIZE = 2000;

    public static double [][] positions = new double[NUM_PARTICLES][2];
    public static double [] forces = new double[NUM_PARTICLES];
    public static double [] densities = new double[NUM_PARTICLES];


    public KdTree neighbours;
    public static List<KdTree.Node> particleCoordinates = new ArrayList<>();
    public static KdTree tree;




    // ------------------------------------------------- //



    @Override
    public void start(Stage stage) throws IOException {

        CONSTANT_FORCE[0] = 0.0;
        CONSTANT_FORCE[1] = -0.1;

        FIGURE_SIZE[0] = 4;
        FIGURE_SIZE[1] = 6;

        Scene scene = new Scene(createContent());

        stage.setTitle("Fluid Dynamics");



        stage.setScene(scene);
        stage.show();


    }

    private Parent createContent(){
        Pane root = new Pane();
        root.setPrefSize(SCREEN_WIDTH,SCREEN_HEIGHT);
        root.setStyle("-fx-background-color: black;");

        for (int i = 0; i < densities.length; i++) {
            densities[i] = 0;
        }

        for (int i = 0; i < NUM_PARTICLES; i++) {
            var p = new Particle();
            particles.add(p);
            root.getChildren().add(p);
        }
        for (int i = 0; i < NUM_PARTICLES/50; i++) {
            for (int j = 0; j < 50; j++) {
                Particle particle = particles.get(i*10+j);
                particle.setTranslateX(i*4+PARTICLE_RADIUS*2+20);
                particle.setTranslateY(j*4+PARTICLE_RADIUS*2+10);

                positions[i*50+j][0] = i*4+PARTICLE_RADIUS*2+20;
                positions[i*50+j][1] = i*4+PARTICLE_RADIUS*2+10;


            }
        }

        for (int i = 0; i < positions.length; i++) {
            particleCoordinates.add(new KdTree.Node(new double[] {positions[i][0],positions[i][1]}));
        }

        tree = new KdTree(2,particleCoordinates);


        /*System.out.println(neighbours.size());
        for (Map.Entry<Double, ArrayList<int []>> entry : neighbours.entrySet()) {
            for (int i = 0; i < entry.getValue().size(); i++) {
                System.out.println("Key: " + entry.getKey());
                for (int j = 0; j < entry.getValue().get(i).length ; j++) {
                    System.out.println("Id " + j + ": " + entry.getValue().get(i)[j]);
                }
            }

        }

         */

        AnimationTimer timer = new AnimationTimer() {
            @Override
            public void handle(long now) {
                onUpdate();
            }
        };
        timer.start();

        return root;
    }

    private static class Particle extends Parent {
        Particle(){
            var shape = new Circle(PARTICLE_RADIUS*2);
            shape.setFill(Color.LIGHTSKYBLUE);

            getChildren().addAll(shape);
        }
    }

    private void onUpdate(){

        KdTree.Node nearest;

        /*for (int i = 0; i < particles.size(); i++) {



            if (ADD_PARTICLES_EVERY == 1){
                if (positions[i][1]+vectorsY[i] > SCREEN_HEIGHT-2*PARTICLE_RADIUS){
                    particles.get(i).setTranslateY(SCREEN_HEIGHT-2*PARTICLE_RADIUS);
                }
                else{
                    particles.get(i).setTranslateY(positions[i][1]+vectorsY[i]);
                }

                particles.get(i).setTranslateX(positions[i][0]+vectorsX[i]);
            }



            if (positions[i][1]+vectorsY[i] > SCREEN_HEIGHT-2*PARTICLE_RADIUS){
                positions[i][1] = SCREEN_HEIGHT-2*PARTICLE_RADIUS;
            }
            else{
                positions[i][1] = positions[i][1]+vectorsY[i];
            }

            positions[i][0] = positions[i][0]+vectorsX[i];

            double currParticle = positions[i][0];
            double currParticleY = positions[i][1];
            nearest = new KdTree.Node(currParticle,currParticleY);
            ArrayList collisions = tree.findCollisions(nearest,PARTICLE_RADIUS);

            System.out.println(collisions);

            //KdTree.Node nearestParticle = tree.findNearest(nearest);
            //System.out.println("Particle at: "+ "("+currParticle+", "+currParticleY+")"+ " is nearest to the particle with id " + "" + " at" +nearestParticle);


            for (int j = 0; j < particles.size(); j++) {
                double otherParticle = positions[j][0];
                double otherParticleY = positions[j][1];
                if (i != j){
                    if (currParticleY >= otherParticleY+PARTICLE_RADIUS & currParticle <=  otherParticle+PARTICLE_RADIUS){

                        if (currParticle <=  otherParticle+PARTICLE_RADIUS){
                            vectorsY[i] = vectorsY[i]-(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                            vectorsY[j] = vectorsY[j]+(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                            vectorsX[i] = vectorsX[i]-(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                            vectorsX[j] =vectorsX[j]+(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                        }
                        else if (currParticle >=  otherParticle+PARTICLE_RADIUS){
                            vectorsY[i] = vectorsY[i]+(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                            vectorsY[j] = vectorsY[j]-(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                            vectorsX[i] = vectorsX[i]+(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                            vectorsX[j] = vectorsX[j]-(PARTICLE_WEIGHT-MOMENTUM_LOSS);
                        }
                    }

                }
            }


            if (currParticle <= PARTICLE_RADIUS){
                vectorsX[i] = -(vectorsX[i]-MOMENTUM_LOSS);
            }
            else if (currParticle >= SCREEN_WIDTH-PARTICLE_RADIUS){
                vectorsX[i] = -(vectorsX[i]-MOMENTUM_LOSS);
            }






            if (currParticleY >= SCREEN_HEIGHT-PARTICLE_RADIUS*2){
                vectorsY[i] = -(vectorsY[i]-MOMENTUM_LOSS)/1.5;
            }

            if (vectorsY[i] > 0){
                vectorsY[i]+=GRAVITATIONAL_FORCE;
            }
            else if (vectorsY[i] < 0){
                vectorsY[i]+=GRAVITATIONAL_FORCE;
            }
            //System.out.println("Particle X: "+ currParticle +" Particle Y:"+currParticleY +" Vector:" + vectorsY[i]);
        }

         */
        for (int i = 0; i < particles.size(); i++) {
            KdTree.Node currParticle = particleCoordinates.get(i);
            ArrayList collisions = tree.findCollisions(currParticle, PARTICLE_RADIUS);
            for (int j = 0; j < collisions.size(); j++) {
                double[] coll_sub = (double[]) collisions.get(j);
                System.out.println(coll_sub[0]);
                densities[j] += NORMALIZATION_DENSITY * (
                        Math.pow(SMOOTHING_LENGTH,2) - getDistance(currParticle.coords_[0],currParticle.coords_[1],coll_sub[0],coll_sub[1])
                        );
            }
            System.out.println();
        }

        //particleCoordinates = new ArrayList<>();
        for (int i = 0; i < positions.length; i++) {
            particleCoordinates.set(i, new KdTree.Node(new double[] {positions[i][0],positions[i][1]}));
        }

        tree = new KdTree(2,particleCoordinates);

        ADD_PARTICLES_EVERY+=1;
        if (ADD_PARTICLES_EVERY == 2){
            ADD_PARTICLES_EVERY = 0;
        }
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
