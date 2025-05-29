import numpy as np
from numpy import pi, sin, cos
import matplotlib.pyplot as plt


class Stick:
    def __init__(self, length, mass, angle, angular_velocity, setpoint):
        self.length = length  # Length of the stick
        self.mass = mass  # Mass of the stick
        self.angle = angle  # Angle of the stick with respect to the vertical (in radians)
        self.angular_velocity = angular_velocity  # Angular velocity of the stick (in radians per second)
        self.gravity = 9.81  # Acceleration due to gravity (m/s^2)
        self.isDampened = True  # Flag to indicate if the stick is dampened


        # PID controller parameters
        self.PID_setpoint = setpoint  # Desired angle (in radians)
        self.previous_error = 0 # Previous error for derivative calculation

        # PID tuning parameters (Ziegler-Nichols method)
        Ku = 100  # Ultimate gain (tuning parameter)         10      5
        Tu = .4 # Ultimate period (tuning parameter)       1.5     11.77

        self.Kp = 0.6*Ku  # Proportional gain
        
        Ti = 0.5*Tu  # Integral time constant (tuning parameter)
        self.Ki = 1.2*Ku/Tu  # Integral gain
        
        Td = 0.125*Tu  # Derivative time constant (tuning parameter)
        self.Kd = 0.075*Ku*Tu # Derivative gain

        if TUNING:
            self.Ki = 0
            self.Kd = 0
            self.Kp = 100
            print(f"PID Tuning Parameters: Kp={self.Kp}, Ki={self.Ki}, Kd={self.Kd}")





    def get_position(self):
        """Calculate the position of the center of mass of the stick."""
        x = self.length / 2 * sin(self.angle)
        y = -self.length / 2 * cos(self.angle)
        return x, y

    def get_torque(self):
        """Calculate the torque acting on the stick."""
        torque = -self.mass * self.gravity * (self.length / 2) * sin(self.angle)
        return torque

    def update_angle(self, dt):
        """Update the angle and angular velocity of the stick."""
        torque = self.get_torque()
        torque += self.PID(dt)  # Add PID control output to the torque
        # Calculate angular acceleration using the torque
        angular_acceleration = torque / (self.mass * (self.length ** 2) / 3)  # Moment of inertia for a rod about its end
        self.angular_velocity += angular_acceleration * dt
        self.angle += self.angular_velocity * dt
        # Keep the angle within the range [0, 2*pi]
        self.angle = self.angle % (2 * pi)

        # Damping to simulate energy loss (optional)
        if self.isDampened:
            damping_factor = 0.99
            self.angular_velocity *= damping_factor
    
    def PID(self, dt):
        """PID controller to control the angle of the stick."""
        integral = 0
        derivative = 0
        
        
        
        error = self.PID_setpoint - self.angle
        
        integral += error * dt
        derivative = (error - self.previous_error) / dt
        self.previous_error = error

        # Proportional term
        p_term = self.Kp * error
        # Integral term
        i_term = self.Ki * integral
        # Derivative term
        d_term = self.Kd * derivative

        # Total PID output
        output = p_term + i_term + d_term
        return output

    
TUNING = False  # Flag to enable or disable PID tuning

# Currently tuned for a 1m stick with a mass of 1kg around the vertical

steps = 2000
dt = 0.01  # Time step (in seconds)
length = 1.0  # Length of the stick (in meters)
mass = 1.0  # Mass of the stick (in kg)
angle = pi - 0.01  # Initial angle (in radians)
angular_velocity = 0.0  # Initial angular velocity (in radians per second)
setpoint =  pi * 11/10   # Desired angle (in radians)

stick = Stick(length, mass, angle, angular_velocity, setpoint)

positions = []
angles = []

for _ in range(steps):
    stick.update_angle(dt)
    positions.append(stick.get_position())
    angles.append(stick.angle)




# Plotting the results
plt.figure(figsize=(10, 5))
plt.subplot(1, 2, 1)
plt.plot([0, stick.length * sin(angle)], [0, -stick.length * cos(angle)], 'b-', label='Initial Position')
#plt.plot([0, stick.length * sin(stick.angle)], [0, -stick.length * cos(stick.angle)], 'r-', label='Final Position', fontsize=18)
plt.plot([0, stick.length * sin(stick.angle)], [0, -stick.length * cos(stick.angle)], 'r-', label='Final Position' )
plt.plot([0, stick.length * sin(stick.PID_setpoint)], [0, -stick.length * cos(stick.PID_setpoint)], 'g--', label='Desired Position' )
plt.plot(0, 0, 'ko', label='Pivot Point')
plt.legend(fontsize=18)
plt.xlim(-1.5, 1.5)
plt.ylim(-1.5, 1.5)
plt.title('Stick Final Position', fontsize=18)
plt.xlabel('X Position (m)', fontsize=18)
plt.ylabel('Y Position (m)', fontsize=18)

plt.grid()
plt.subplot(1, 2, 2)

plt.plot(np.arange(0, steps * dt, dt), angles, 'b.')
plt.axhline(y=stick.PID_setpoint, color='g', linestyle='--', label='Desired Angle')
plt.axhline(y=0, color='k', linestyle='--', label='Vertical')
plt.axhline(y=angle, color='b', linestyle='--', label='Initial Position')
plt.legend(fontsize=18)
plt.title('Angle vs Time', fontsize=18)
plt.xlabel('Time (s)', fontsize=18)
plt.ylabel('Angle (radians)', fontsize=18)
plt.grid()
plt.tight_layout()



if TUNING:
    plt.subplot(1, 3, 3)
    plt.title('PID Closeup Tuning')
    plt.plot(np.arange(0, steps * dt, dt), angles, 'b.')
    plt.axhline(y=stick.PID_setpoint, color='g', linestyle='--', label='Desired Angle')
    plt.legend()
    plt.xlabel('Time (s)')
    plt.ylabel('Angle (radians)')
    plt.grid()
    plt.ylim(0.98*setpoint, 1.02*setpoint)
 
plt.show()


