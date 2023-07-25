package com.example.fluiddynamicsjavafx;

import java.io.BufferedReader;
import java.io.IOException;
import java.io.InputStreamReader;
import java.io.PrintWriter;
import java.net.Socket;
import java.util.ArrayList;
import java.util.List;

public class Client implements Runnable{

    private Socket client;
    private BufferedReader in;
    private PrintWriter out;
    private boolean done;
    List<Particle> particles;

    public double[] get_doubleArr(String p){
        p = p.replaceAll("[()]", "");
        p = p.replace("[", "");
        p = p.replace("]", "");
        String [] pArr = p.split(",");
        double[] result = new double[pArr.length];
        for (int i = 0; i < pArr.length; i++) {
            result[i] = Double.parseDouble(pArr[i]);
        }
        return result;
    }

    public void set_particles(String[] data){
        boolean empty = false;
        if (particles.size() == 0){
            System.out.println("-----------------------------------------------------------------");
            System.out.println("It's empty!");
            System.out.println("-----------------------------------------------------------------");
            empty = true;
        }

        for (String s : data){
            if (!s.equals("/update")){
                //System.out.println(s);
                String[] point = s.split("_");
                double[] coords = get_doubleArr(point[0]);
                int id = Integer.parseInt(point[1]);
                double density = Double.parseDouble(point[2]);
                double pressure = Double.parseDouble(point[3]);
                double[] forces = get_doubleArr(point[4]);
                if (empty){
                    Particle p = new Particle(coords,id);
                    p.setDensity(density);
                    p.setPressure(pressure);
                    p.setVelocity(forces);
                    particles.add(p);
                }
                else {
                    Particle p = particles.get(id);
                    p.setDensity(density);
                    p.setPressure(pressure);
                    p.setVelocity(forces);
                    p.setCoords_(coords);
                }

            }
        }
    }

    @Override
    public void run() {
        try {
            client = new Socket(Fluid2D.hostAddress, Fluid2D.hostPort);
            out = new PrintWriter(client.getOutputStream(),true);
            in = new BufferedReader(new InputStreamReader(client.getInputStream()));
            particles = new ArrayList<>();
            InputHandler inputHandler = new InputHandler();
            Thread t = new Thread(inputHandler);
            t.start();

            String inMessage;
            while ((inMessage = in.readLine()) != null){
                String[] data = inMessage.split(" ");
                if (inMessage.startsWith("/calculate")){
                    // extract data and calculate, send back calculated data
                    ArrayList<Particle> points = new ArrayList<>();
                    int i = Integer.valueOf(data[1]);
                    int num_connections = Integer.valueOf(data[2]);
                    for (int j = 0; j < particles.size(); j++) {
                        if ((particles.get(j).getCoords_()[0] > (Fluid2D.MAX_X*i/num_connections) || (i==0 && particles.get(j).getCoords_()[0] >= ((Fluid2D.MAX_X*i/num_connections))) ) && particles.get(j).getCoords_()[0] <= (Fluid2D.MAX_X*(i+1)/num_connections)){
                            points.add(particles.get(j));
                        }
                    }
                    out.println("/wd "+ calculateSegment(points));

                }
                else if (inMessage.startsWith("/update ")){
                    // update the particles
                    set_particles(data);
                }
                else {
                    System.out.println("Unrecognized command. "+ inMessage);
                }
            }

        }catch (IOException e){
            //handle here
        }
    }

    public void shutdown(){
        done = true;
        try{
            in.close();
            out.close();
            if (!client.isClosed()){
                client.close();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }

    }



    class InputHandler implements Runnable{

        @Override
        public void run() {
            try {
                BufferedReader inputReader = new BufferedReader(new InputStreamReader(System.in));
                while (!done){

                    String message = inputReader.readLine();
                    String[] data = message.split(" ");
                    if (message.equals("/quit")){
                        inputReader.close();
                        shutdown();
                    }
                    else {
                        out.println(message);
                    }
                }
            }catch (IOException e){
                //handle here
            }
        }
    }

    public List<Particle> calculateSegment(List<Particle> particles){
        if (Fluid2D.grid == null){
            Fluid2D.grid = new Grid(Fluid2D.size,Fluid2D.smoothingLength);
            for (Particle p : particles){
                Fluid2D.grid.Insert(p);
            }
        }

        for(Particle particle : particles){
            Fluid2D.neighbors[particle.id] = (Fluid2D.grid.GetNeighbors(particle));
        }

        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            Fluid2D.calculateDensity(currParticle);
            if (currParticle.getDensity() == 0){
                currParticle.setDensity(0.0001);
            }
        }

        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            Fluid2D.calculatePressure(currParticle);
        }



        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            Fluid2D.calculateForces(currParticle);
        }
        for (int i = 0; i < particles.size(); i++) {
            Particle currParticle = particles.get(i);
            currParticle.setVelocity(Fluid2D.newForces[currParticle.id]);
        }


        for (Particle currParticle : particles){
            currParticle.setVelocity(Fluid2D.newForces[currParticle.id]);
            currParticle.setCoords_(Fluid2D.add(currParticle.getCoords_(), Fluid2D.mul(currParticle.getVelocity(),Fluid2D.timeStep)));
        }
        for (int i = 0; i < particles.size(); i++) {
            checkIfBounced(particles.get(i));
        }



        return particles;
    }

    public static void checkIfBounced(Particle currparticle){
        double []coords = currparticle.getCoords_();
        double []velocity = currparticle.getVelocity();
        if (coords[0] <= Fluid2D.MIN_X){

            double [] newVelocity = {-velocity[0]*Fluid2D.COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {0, currparticle.getCoords_()[1]});
        }
        if (coords[0] >= Fluid2D.MAX_X){

            double [] newVelocity = {-velocity[0]*Fluid2D.COR,velocity[1]};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] {Fluid2D.MAX_X, currparticle.getCoords_()[1]});
        }
        if (coords[1] <= 0){
            //double [] newVelocity = {velocity[0],-velocity[1]};
            double [] newVelocity = {velocity[0],-velocity[1]*Fluid2D.COR};
            currparticle.setVelocity(newVelocity);
            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],0});
        }
        if (coords[1] >= Fluid2D.MAX_Y){

            double [] newVelocity = {velocity[0],-velocity[1]*Fluid2D.COR};
            //double [] newVelocity = {velocity[0],-velocity[1]};
            currparticle.setVelocity(newVelocity);

            currparticle.setCoords_(new double[] { currparticle.getCoords_()[0],Fluid2D.MAX_Y});
        }}

    public static void main(String[] args) {
        Thread clientThread = new Thread(new Client());
        clientThread.start();
    }

}
