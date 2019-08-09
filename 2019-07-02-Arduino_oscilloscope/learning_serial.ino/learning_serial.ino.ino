int analogPin0 = A0; // potentiometer wiper (middle terminal) connected to analog pin 3
int analogPin1 = A1;

                    // outside leads to ground and +5V
int val0 = 0;  // variable to store the value read
int val1 = 0;  // variable to store the value read
int val2 = 0;

void setup() {
  Serial.begin(115200);           //  setup serial
}

void loop() {
  //Square Wave
  for(int a=0; a<=100; a++)
  {
    val0 = analogRead(analogPin0);  // read the input pin
    val1 = analogRead(analogPin1);
    val2 = 0;
    Serial.print(val0);        
    Serial.print(" ");
    Serial.print(val1);
    Serial.print(" ");
    Serial.println(val2);
    delay(1);  
  }
  for(int a=0; a<=100; a++)
  {
    val0 = analogRead(analogPin0);  // read the input pin
    val1 = analogRead(analogPin1);
    val2 = 100;
    Serial.print(val0);         
    Serial.print(" ");
    Serial.print(val1);
    Serial.print(" ");
    Serial.println(val2);
    delay(1);
  }
}
