# Generated by Django 3.0.5 on 2020-04-03 13:08

from django.db import migrations


class Migration(migrations.Migration):

    dependencies = [
        ('drugstone', '0020_auto_20200403_1302'),
    ]

    operations = [
        migrations.RemoveField(
            model_name='protein',
            name='proteins',
        ),
    ]
