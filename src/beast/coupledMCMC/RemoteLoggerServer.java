package beast.coupledMCMC;

import java.io.DataInputStream;
import java.io.EOFException;
import java.io.File;
import java.io.FileOutputStream;
import java.io.IOException;
import java.io.PrintStream;
import java.net.ServerSocket;
import java.net.Socket;

import beast.core.util.Log;

public class RemoteLoggerServer {

	// server is listening on port
	private ServerSocket ss;
	private Thread serverThread;

	public RemoteLoggerServer(int port) throws IOException {
		// server is listening on port
		ss = new ServerSocket(port);

		serverThread = new Thread() {
			boolean running = true;

			@Override
			public void run() {
				// running infinite loop for getting client requests
				while (running) {
					Socket s = null;

					try {
						// socket object to receive incoming client requests
						s = ss.accept();

						Log.info("A new client is connected : " + s);

						// obtaining input and out streams
						DataInputStream dis = new DataInputStream(s.getInputStream());

						// create a new thread object
						Thread t = new LoggerClientHandler(s, dis);

						// Invoking the start() method
						t.start();

					} catch (Exception e) {
						try {
							s.close();
						} catch (IOException e1) {
							// TODO Auto-generated catch block
							e1.printStackTrace();
						}
						e.printStackTrace();
					}
				}
			}
			
			@Override
			public void interrupt() {
				running = false;
				super.interrupt();
			}
		};
			
		serverThread.start();

	}

	public class LoggerClientHandler extends Thread {
		final DataInputStream dis;
		final Socket s;
		PrintStream out = null;

		// Constructor
		public LoggerClientHandler(Socket s, DataInputStream dis) {
			this.s = s;
			this.dis = dis;
		}

		@Override
		public void run() {
			String received;

			// first message indicates the file name
			try {
				received = dis.readUTF();
				Log.warning("Creating file at " + received);
				if (received.equals("null")) {
					out = System.out;
				} else {
					if (new File(received).exists()) {
						out = new PrintStream(new FileOutputStream(received, true));
					} else {
						out = new PrintStream(received);
					}
				}
			} catch (IOException e1) {
				// TODO Auto-generated catch block
				e1.printStackTrace();
			}

			while (true) {
				try {

					// receive the answer from client
					received = dis.readUTF();
					if (received.equals("close")) {
						break;
					}
					if (received.length() > 0) {
						out.println(received);
						out.flush();
					}
				} catch (EOFException e) {
					break;
				} catch (IOException e) {
					e.printStackTrace();
				}
			}

			try {
				// closing resources
				this.dis.close();
				// close connections
				serverThread.interrupt();
				ss.close();
			} catch (IOException e) {
				e.printStackTrace();
			}
		}
	}
}
