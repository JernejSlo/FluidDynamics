package com.example.fluiddynamicsjavafx;

import java.io.*;
import java.net.ServerSocket;
import java.net.Socket;
import java.util.ArrayList;
import java.util.List;
import java.util.concurrent.ExecutorService;
import java.util.concurrent.Executors;

public class Server implements Runnable{

    public ServerSocket server;
    private ArrayList<ConnectionHandler> connections;
    private boolean done;
    private ExecutorService pool;
    private ArrayList<Boolean> calculated;
    public List<Particle> particleCoordinates;

    public Server(List<Particle> particleCoordinates){
        connections = new ArrayList<>();
        done = false;
        this.particleCoordinates = particleCoordinates;
    }

    @Override
    public void run() {
        try {
            ServerSocket serverSocket = new ServerSocket(Fluid2D.hostPort);
            server = serverSocket;
            pool = Executors.newCachedThreadPool();
            calculated = new ArrayList<>();
            while (!done){
                System.out.println(server);
                Socket client = server.accept();
                ConnectionHandler connectionHandler = new ConnectionHandler(client);
                connections.add(connectionHandler);
                pool.execute(connectionHandler);
            }
        } catch (IOException e) {
            shutdown();
        }
    }

    public void broadcast(String message){
        for (ConnectionHandler ch : connections){
            if (ch != null){
                ch.sendMessage(message);
            }
        }
    }

    public void broadcastI(String message, String additional){
        for (ConnectionHandler ch : connections){
            if (ch != null){
                ch.sendMessage(message + " " + ch.id + " " + additional);
            }
        }
    }

    public void shutdown() {
        try {
            done = true;
            if (server != null && !server.isClosed()) {
                server.close();
            }
            for (ConnectionHandler ch : connections) {
                ch.shutdown();
            }
        } catch (IOException e) {
            e.printStackTrace();
        }
    }


    class ConnectionHandler implements Runnable{

        private Socket client;
        private BufferedReader in;
        private PrintWriter out;
        private String message;


        private int segment;
        public int id;

        public ConnectionHandler(Socket client){
            this.client = client;
        }

        public boolean done(){
            boolean allTrue = true;
            for (boolean value : calculated) {
                if (!value) {
                    allTrue = false;
                    break;
                }
            }
            return allTrue;
        }

        public void reset(){
            for (int i = 0; i < calculated.size(); i++) {
                calculated.set(i,false);
            }
        }

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

            for (String s : data){
                if (!s.equals("/wd")){
                    String[] point = s.split("_");
                    double[] coords = get_doubleArr(point[0]);
                    int id = Integer.parseInt(point[1]);
                    double density = Double.parseDouble(point[2]);
                    double pressure = Double.parseDouble(point[3]);
                    double[] forces = get_doubleArr(point[4]);

                    Particle p = particleCoordinates.get(id);
                    p.setDensity(density);
                    p.setPressure(pressure);
                    p.setVelocity(forces);
                    p.setCoords_(coords);
                }
            }
        }

        @Override
        public void run() {
            try {
                out = new PrintWriter(client.getOutputStream(), true);
                in = new BufferedReader(new InputStreamReader(client.getInputStream()));
                out.println("Lahko napišeš komando.");
                id = connections.size()-1;
                calculated.add(true);
                if (calculated.size() >= 1){
                    broadcast("/update " + particleCoordinates);
                    broadcastI("/calculate",""+calculated.size());
                }


                while((message = in.readLine()) != null){
                    if (message.startsWith("/wd")){
                        //System.out.println("data from "+ id +" has been sent from client to server!");
                        calculated.set(id, true);
                        String[] data = message.split(" ");
                        if (done()){
                            set_particles(data);
                            reset();
                            Fluid2D.updateParticlesDistributed();
                            broadcast("/update " + particleCoordinates);
                            broadcastI("/calculate",""+calculated.size());
                        }
                        // zapiši podatke
                        // preveri če je calculated cel končan in ga spremeni na false
                        // pošlji komando da zračunajo naslednji del
                    }
                    else if (message.startsWith("/disconnect")){
                        System.out.println("A client has disconnected");
                        String[] data = message.split(" ");
                        calculated.remove(Integer.valueOf(data[1]));
                        if ( Integer.valueOf(data[1]) < id){
                            id-= 1;
                        }
                    }
                    else{
                        //HUH? ni komande pac
                        System.out.println("Client sent invalid message. " + message);
                    }
                }
            }
            catch (IOException e){
                shutdown();
            }
        }

        public void sendMessage(String message){
            out.println(message);
        }

        public void shutdown(){

            if (!client.isClosed()){
                try {
                    in.close();
                    out.close();
                    client.close();
                } catch (IOException e) {
                    e.printStackTrace();
                }
            }
        }

    }

}
